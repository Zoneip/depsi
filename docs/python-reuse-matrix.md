# 三方 Python 代码差异与复用清单

状态：**分析清单（2026-10-29）**，本轮不改代码。范围：`depsipy`（当前实现）、
`reference/python`（官方 TUDelftGeodesy/DePSI）、`workspace/Doris干涉加时序depsi`（RapidSAR/DORIS 干涉链）。
许可边界（**2026-10-29 用户更新**）：官方 Python（Apache-2.0）与 **Doris 干涉链代码均可直接使用**；
Doris 侧原先“许可不明、只参考”的结论作废。

## 0. 结论摘要

1. `depsipy` 继续作为主架构（numpy + dataclass + raw memmap），不引入 xarray/Dask。
2. 官方 Python 版是 **xarray/Dask + 多依赖**架构，整体不能直接合并；只适配其中
   **自包含、纯 numpy/scipy** 的函数与算法。
3. 四个新功能中：
   - **F1 稠密化**：官方版只支持 `n_connections=1`，缺 `Nest/2of3/3of3`，只能复用它
     的 KDTree 连接与模糊度合成思路，主体按 MATLAB + `SPEC-003` 实现。
   - **F2 大气滤波+Kriging**：官方版算法完整，但强依赖 xarray/pykrige；适配滤波权重、
     变差函数与 Kriging 方程到 numpy/sparse，不直接调用官方函数。
   - **F3 噪声滤波+STC**：官方版有 STC 公式与 `find_points_within_buffer` 可适配；
     噪声滤波官方版没有，按 MATLAB `ps_noise_filter` 重写。
   - **F4 输出/地理编码**：官方 `transformations.py` 是纯 numpy/scipy，可复用价值最高；
     Doris 的 geocoding 只作流程参考。
4. Doris 链不含时序反演，F1/F2/F3 基本不相关；其价值在 F4 的裁剪/重采样/地理编码流程参考。

## 1. 三方源码盘点

| 来源 | 路径 | 规模 | 角色 | 许可 | 主要依赖 |
| --- | --- | --- | --- | --- | --- |
| depsipy | `src/depsipy` | 14 文件 / 3,912 行 | 当前主实现 | DePSI 科研许可 | numpy、scipy、numba（可选） |
| 官方 Python | `reference/python/depsi` | 21 文件 / 11,854 行 | 官方 beta 重实现 | **Apache-2.0** | xarray、Dask、pykrige、sarxarray、stmtools、ruptures、networkx、drama、asf_search、scikit-learn、sparse |
| Doris 链 | `workspace/Doris干涉加时序depsi/doris干涉` | 150 文件 / 79,241 行 | 干涉链（无时序 PS 反演） | 用户 2026-10-29 确认**可直接使用** | DORIS 可执行文件、外部 SAR 数据、Python 脚本 |

官方测试：`reference/python/tests/` 20 个测试文件，覆盖 densification、atmosphere、
point_quality、transformations 等，可作为等价性参考（Apache-2.0）。

## 2. 架构与数据模型差异

| 维度 | depsipy | 官方 Python | Doris |
| --- | --- | --- | --- |
| 数据模型 | numpy 数组 + dataclass；像素坐标 1-based | `xr.Dataset/DataArray`，`space/time` 维、`pnt_uid`、attrs | 文件驱动；DORIS 中间产物 |
| I/O | raw（memmap、无头、复像素交错）、`npz/json` 状态 | HDF5/Zarr/NetCDF 风格，`stm` 数据集 | DORIS 二进制/文本 + GeoTIFF |
| 执行 | 单进程 numpy/numba；函数级 | Dask 并行、chunking | shell/Python 批处理，调外部程序 |
| 坐标 | az/r 像素 + spacing | `azimuth/range` 或 `x/y` 坐标，单位在 attrs | 雷达坐标 + 地理坐标 |
| 参数 | MATLAB `param_file.txt` | pydantic 配置（`config.py`） | 脚本内/配置文件 |
| 依赖重量 | 轻 | 重（xarray/Dask/pykrige 等） | 重（DORIS 环境） |

因此：**官方代码不能整包 vendor**；可复用函数必须去掉 xarray 接口，落到 numpy 数组。

## 3. 四个功能的差异与复用矩阵

图例：✔ 已有；◐ 部分；✘ 未实现；↺ 可适配复用；◎ 仅参考。

### F1 · step6 稠密化（需求 R1 首要项）

| 来源 | 现状 | 可复用内容 | 差距 |
| --- | --- | --- | --- |
| MATLAB | `ps_densification.m` 1282 行，21 for + 1 while；`Nest=3`、`2of3/3of3/mode_fix`、幅度分箱迭代 | 算法语义与验收口径 | 尚未实现 |
| 官方 Python | [densification.py](file:///home/zoneip/工程/Zoneip/depsi/reference/python/depsi/densification.py#L11-L221) `densification()`、`_determine_dens_connections()` | KDTree 最近邻连接、`h2ph/dd_phase` 构造、弧模糊度→点模糊度合成 | **只支持 `n_connections=1`**；无 `Nest/2of3/3of3`；用 periodogram，不是 ILS；无幅度分箱/迭代 |
| Doris | 无时序反演 | — | 不相关 |
| depsipy | ✘；已有 [SPEC-003](file:///home/zoneip/工程/Zoneip/depsi/specs/SPEC-003-稠密化向量化设计.md) 设计 | `estimate_arcs` 批量 ILS、`form_network`/`unwrap_network` | 需实现 |

**复用决定**：↺ 只适配 `_determine_dens_connections` 的 KDTree 思路与“网络点模糊度 + 弧模糊度”
合成公式；`2of3/3of3/mode_fix`、`Nest`、幅度分箱按 MATLAB/SPEC-003 实现。
**测试**：官方 `tests/test_densification.py` 可提供 KDTree/连接用例；主体用 SPEC-003 的
小数据对拍（PS 数、acheck、参数、ens_coh）。

### F2 · step4 大气滤波 + Kriging

| 来源 | 现状 | 可复用内容 | 差距 |
| --- | --- | --- | --- |
| MATLAB | `ps_filtering.m`、`ps_atmo_kriging.m`、`ps_variogram`/`ps_fit_vario`/`ps_kriging` | 语义、参数、验收 | 未实现 |
| 官方 Python | [atmosphere_estimation.py](file:///home/zoneip/工程/Zoneip/depsi/reference/python/depsi/atmosphere_estimation.py#L59-L150) `estimate_unmodeled_displacement`；L153-255 变差函数；L372-478 `setup_kriging_system`；L481-725 `solve_kriging`；L781+ `estimate_atmosphere_phase` | 高斯/三角/block 滤波权重矩阵；经验变差函数分箱/拟合；Kriging 方程与预测 | 强依赖 xarray/Dask/pykrige；MATLAB 用自定义 variogram + 普通 Kriging，官方用 universal Kriging |
| Doris | 干涉链滤波/去平地 | 流程参考 | 无时序大气分离 |
| depsipy | ✘ | 无 | 需新增 `steps/filtering.py`、`steps/kriging.py` |

**复用决定**：↺ 复用滤波权重与变差函数数学（改写为 numpy，输入 `(Nps, Nifgs)` 数组）；
`setup_kriging_system`/`solve_kriging` 只作算法参考，用 scipy/sparse 实现普通 Kriging，
**不新增 pykrige 依赖**。
**测试**：官方 `tests/test_atmosphere_estimation.py` 的滤波/变差/Kriging 用例；
另加合成 APS（已知空间相关长度/幅度）验收。

### F3 · step8/9 噪声滤波 + STC

| 来源 | 现状 | 可复用内容 | 差距 |
| --- | --- | --- | --- |
| MATLAB | `ps_noise_filter.m`（7 for）、`ps_spatio_temporal_consistency.m`（8 for + 1 parfor） | 语义、`stc_min_max`、`filter_length` | 未实现 |
| 官方 Python | [point_quality.py](file:///home/zoneip/工程/Zoneip/depsi/reference/python/depsi/point_quality.py#L805-L866) `compute_spatiotemporal_consistency`；`find_points_within_buffer`（network.py） | STC 公式（时间差分标准差、min、mm 换算）、距离缓冲查询 | 逐点 Python for；官方**没有** `ps_noise_filter` 等价函数 |
| Doris | 无 | — | 不相关 |
| depsipy | ✘ | `sidelobe_mask`/selection 过滤模式可参考 | 需新增 |

**复用决定**：↺ 适配 STC 公式与 buffer 查询（用 `cKDTree`/广播向量化）；噪声滤波按 MATLAB
`ps_noise_filter` 重写，官方版无对应实现。
**测试**：官方 `tests/test_point_quality.py` 的 STC 用例；新增合成时序滤波/STC 验收。

### F4 · step7 输出 + 地理编码

| 来源 | 现状 | 可复用内容 | 差距 |
| --- | --- | --- | --- |
| MATLAB | `ps_output.m`/`ps_output_nodens.m`；`geocoding_v0.8` | 产物格式、resfile、地理编码流程 | 未实现 |
| 官方 Python | [transformations.py](file:///home/zoneip/工程/Zoneip/depsi/reference/python/depsi/transformations.py#L46-L703) `latlonh_to_xyz`/`xyz_to_latlonh`/`orbit_fit`/`xyz_to_topocentric`/`radar_to_latlonh`/`latlonh_to_radar` 等 | **纯 numpy/scipy 坐标/轨道变换**，强复用候选；`visualization.py` 仅绘图 | 官方**没有** raw/GeoTIFF 产物写出；数据模型是 xarray |
| Doris | `RapidSAR/Script/Geocoding_*.py`、`sar_image_crop.py`、`*resample*` | MFF 头/`.res`/`lph2xyz` 流程已移植 | 依赖 DORIS 可执行文件与外部数据；代码可直接使用 |
| depsipy | ✘；已有 raw 读写 | `io/raw.write_raw`、`steps/unwrapping` 结果 | 需新增 `steps/output.py` / `steps/geocoding.py` |

**复用决定**：↺ 适配 `transformations.py` 的坐标/轨道函数（保留 Apache 头与 NOTICE）；
输出/GeoTIFF 由 depsipy 实现（可选 `rasterio`，不进核心依赖）；Doris geocoding 的
`metadataNew/intrp_orbit/lph2xyz` 已按算法移植到 `io/doris.py` + `core/orbit.py` +
`steps/geocoding.py`（无 GDAL 依赖）。
**测试**：官方 `tests/test_transformations.py` 的坐标正反算用例；raw 产物用逐字节 golden 校验。

## 4. 模块级复用清单

| 官方文件/函数 | depsipy 目标 | 动作 | 许可 | 风险 |
| --- | --- | --- | --- | --- |
| `densification._determine_dens_connections` | `steps/densification.py` 邻接辅助 | ↺ 适配 | Apache-2.0 | 只处理 1 连接；需扩展到 `Nest` |
| `densification.densification` | 仅参考 | ◎ | Apache-2.0 | 不满足 `2of3/3of3`；不整函数复制 |
| `atmosphere_estimation.estimate_unmodeled_displacement` | `steps/filtering.py` | ↺ 适配 | Apache-2.0 | 权重矩阵语义需与 MATLAB filter_nl 对齐 |
| `calculate_variogram_cloud` / `calculate_empirical_variogram` / `fit_variogram` | `steps/kriging.py` | ↺ 适配 | Apache-2.0 | 分箱/拟合口径需对拍 |
| `setup_kriging_system` / `solve_kriging*` / `estimate_atmosphere_phase` | 仅参考 | ◎ | Apache-2.0 | pykrige/xarray 依赖不引入 |
| `point_quality.compute_spatiotemporal_consistency` | `steps/stc.py` | ↺ 适配 | Apache-2.0 | 逐点循环需向量化 |
| `find_points_within_buffer`（`network.py`） | `steps/stc.py` 邻域查询 | ↺ 适配 | Apache-2.0 | 地理/欧氏模式 |
| `postprocessing.stm_point_filter` | 通用点过滤辅助 | ↺ 适配 | Apache-2.0 | 极简函数，可直接改写 |
| `transformations.*`（坐标/轨道） | `steps/geocoding.py` | ↺ 适配 | Apache-2.0 | 需要 `master_res`/轨道元数据 |
| `visualization.plot_*` | 可选绘图 | ◎ | Apache-2.0 | 不进入核心链 |
| Doris `Geocoding_*.py` / `sar_image_crop.py` / `metadataNew.m` | `io/doris.py`、`core/orbit.py`、`steps/geocoding.py` | ● 移植 | 可直接使用 | 仅 numpy/scipy，不引入 GDAL/DORIS |

## 5. 许可与依赖策略

1. 直接复用/适配的 Apache-2.0 代码：保留原版权头，新增 `THIRD_PARTY_NOTICES.md`，
   标明来源仓库、commit `51939cb` 与修改说明。
2. Doris/RapidSAR 代码（2026-10-29 用户确认可直接使用）：移植处在文件头注明来源
   （如 `io/doris.py` ← `metadataNew.m`/MFF 约定，`core/orbit.py` ← `intrp_orbit.m`，
   `steps/geocoding.py` ← `lph2xyz.m`/`xyz2t.m`/`xyz2lp.m`）；不引入 GDAL/DORIS 运行时依赖。
3. 不引入 xarray/Dask/pykrige/sarxarray 等重依赖到 depsipy 核心；Kriging 用 scipy/sparse
   实现；GeoTIFF 作为可选 extra（`rasterio`），不默认安装。
4. 官方测试用例若被改写进 depsipy，需在测试文件头注明 Apache-2.0 来源。

## 6. W 型迭代 loop（本轮为分析，后续按功能实现）

| 批次 | 需求 ↔ 验收 | 设计 ↔ 差异核对 | 实现 ↔ 单测 | 集成 ↔ 基准/文档 |
| --- | --- | --- | --- | --- |
| W1（本轮） | 四个功能都有来源/差距/复用结论 | 本清单 + 源码哈希 | 不写代码 | 更新 `docs/README.md` 索引 |
| W2 F1 稠密化 | SPEC-003 验收 | KDTree/模糊度合成与官方对比 | 官方 densification 用例 + MATLAB 小样本 | step6 小数据基准 |
| W3 F2 大气 | APS 残差/相关系数验收 | 滤波/变差/Kriging 公式核对 | 官方 atmosphere 用例 + 合成 APS | step4 基准 |
| W4 F3 噪声/STC | STC 分布/滤波窗口验收 | STC 公式与官方逐条对齐 | 官方 point_quality 用例 | step8 基准 |
| W5 F4 输出/地理编码 | 产物 golden 校验 | 坐标正反算与官方对比 | 官方 transformations 用例 | 输出互通测试 |

终点：四个功能全部有“来源 → 复用决定 → 测试用例 → 基准入口”，再进入实现阶段。
本轮不做代码合并，不引入新依赖，不改公开 API。

## 7. 实施进度（2026-10-29，第一轮）

| 功能 | 已落地 | 测试 | 尚未完成 |
| --- | --- | --- | --- |
| F1 稠密化 | `steps/densification.py` + `pipeline.run_densification` | `tests/test_densification.py`（6 项）+ `test_pipeline.py` | 幅度分箱迭代、`mode_fix`、最终模型合并 |
| F2 大气 | `steps/filtering.py`、`steps/kriging.py`：滤波权重、变差函数、普通 Kriging | `tests/test_filtering_kriging.py`（6 项） | universal Kriging、局部块优化、端到端 |
| F3 噪声/STC | `steps/noise_filter.py`、`steps/stc.py` + pipeline | `tests/test_noise_stc.py`（3 项）+ `test_pipeline.py` | 网格缓冲分块（当前 KD-tree 等价实现） |
| F4 输出/地理编码 | `steps/output.py`、`steps/geocoding.py`、`pipeline.run_output` | `tests/test_output_geocoding.py`（3 项）+ `test_pipeline.py` | `_defo`/`_trend` 整幅栅格、雷达/轨道反演 |

全量测试：124 项（124 passed，无 strict xfail；云端 122 passed + 2 skipped）。

