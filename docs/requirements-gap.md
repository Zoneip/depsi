# 需求 ↔ 实现 ↔ 验证 对照

状态：**现行（2026-10-29）**。需求来源为用户提供的两份材料；行为来源为实际源码。
`待确认` 条目不得当作结论使用。

## 0. 需求原文（`workspace/Depsi效率提升需求.docx`）

> 1. step3、step6（稠密化）等步骤效率很低，是主要的优化步骤，数据增长会导致效率快速降低，
>    无法满足工程需求。**（重要）**
> 2. 对于后置的结果处理尽量调通 8、9 步，可以实现对结果的时间维平滑。
> 3. 建议在小数据测试有明显提升，在扩大研究区数据，以此减少试错成本。

参数说明来自 `workspace/Depsi软件参数.xlsx`（64 行，含单位与取值）；等效 Markdown 清单已
迁移到 [`rules/01-结构.md`](../rules/01-结构.md) 的「工程参数清单」一节。

## 1. 业务步骤 ↔ 源码实际调用

| 用户术语 | 源码位置 | 实际调用的函数 | 备注 |
| --- | --- | --- | --- |
| 第 1 步 | `depsi.m` `processing_groups==1` | 读参数、读输入清单、建立常量 | 也处理 `breakpoint>1000` 时按 orbit 反查索引 |
| 第 2 步 | `==2` | `ps_calibration`、`ps_selection` | 产出 `_psc_2orig_sel*`、`_psp_2orig_sel*`、`_mrm`、`_amp_disp` |
| 第 3 步 | `==3` | `ps_network`（内部 `ps_vce`→`ps_form_network`→`ps_estimation`→`ps_spatial_unwrapping`） | 去趋势后的第一轮组网/估计/解缠 |
| 第 4 步 | `==4` | 再次 `ps_network`、`ps_filtering`、`ps_atmo_kriging` | 大气与形变分离 |
| **稠密化（用户称 step6）** | `==5` **和** `==6` | `ps_densification`（`depsi.m:909` 与 `:1098`） | `==6` 整体被 `if ~isempty(defo_method)` 包住；参数表 `defo_method=[]` ⇒ **该分支默认不执行** |
| 第 7 步 | `==7` | `ps_spatio_temporal_consistency(...,'orig',...)` | 输出与一致性 |
| **8 步 / 9 步（时间维平滑）** | `==8` | `ps_noise_filter` + `ps_spatio_temporal_consistency(...,'filt',...)` | 参考源码**只有 1..8 组**，没有第 9 组 |

结论：`step6` 与 `step9` 是**生产侧的编号口径**，与本仓库参考源码的组号不一致。
生产入口未提供前，一律按**函数名**（`ps_densification`、`ps_noise_filter`、
`ps_spatio_temporal_consistency`）标识任务，不猜编号。

## 2. 功能矩阵（Python 侧 depsipy）

图例：`✔` 已实现且有证据；`◐` 部分实现或未对拍；`✘` 未实现；`?` 待生产入口确认。

| 能力 | MATLAB 函数 | depsipy | 证据 / 缺口编号 |
| --- | --- | --- | --- |
| raw 栅格读写（cpx float32/int16、裁剪、稀疏像素） | `freadbk`/`fwritebk` | ✔ | `tests/test_io_raw.py`（8 项） |
| Doris 产物读取（MFF 头、`.r00/.x00/.j00`、`.res` 元数据/轨道） | Doris/RapidSAR | ✔ | `io/doris.py`；`tests/test_doris_geocoding.py`（真实 `example.res` 解析 + 合成 MFF 往返） |
| 参数文件解析 | `readinput`/`ps_readinput_parameters` | ✔ | `tests/test_io_params.py`（9 项） |
| 输入清单解析与排序 | `depsi.m` 读取段 | ✔ | `tests/test_simulation_selection.py`；预检另做结构校验 |
| 27 个人工模型 + 断点 | `ps_model_definitions` | ✔ | `test_model_definitions_match_matlab_table` 等 |
| LAMBDA/ILS/bootstrap | `decorrel_freek`/`ps_chistart`/`bk_lsearch_count` | ✔ | 与 Octave 逐量一致（`docs/performance.md` §2） |
| 几何/OMT/相干性 | `ps_ils_and_bootstrap_estimation` | ✔ | 同上；`varfac` max|Δ|=1.8e-14 |
| 周期图估计（`ps_method='perio'`） | `ps_periodogram_estimation` | ✘ | 预检报 `python_gap` |
| PSC 网格选择（`threshold`） | `ps_selection` | ◐ | 选点精度 1.00、与 v1.7.8.8 逐个相同；但 `_grid_select` 未含 `eachgrid`/`coherence_eachgrid` 分支 |
| `psc_distribution='uniform'` | `ps_selection` | ◐ | 有实现，**未与 MATLAB 对拍** |
| 旁瓣掩膜（+hhf 去簇） | `ps_selection` | ✔ | 已修 `GAP-S02`：扫描改为“最强为参考 + 连续相干邻域”，并修掉 `sidelobe_mask` 的索引/掩膜误用；同数据 Octave 对拍 PSC/PSP 集合完全一致（`scripts/octave_parity.py`、`test_matlab_parity.py`） |
| PSP 输出载荷 | `ps_selection` | ✔ | 已修 `GAP-S01`：`_psp_2orig_sel*` 写 `[cell_az cell_r az r phase h2ph]`，`_psp_validation_sel*` 写 `[amp_disp amp]`；`test_psp_export_has_observation_payload` |
| 幅度/相干性校准 | `ps_calibration` | ✘ | 仅支持 `calfactors=1` |
| 弧网（Delaunay） | `ps_form_network` | ✔ | 选择/网络测试 |
| spider 冗余网 | `ps_construct_redundant_network` | ◐ | 用 kd-tree 近似，**未复现 Nparts 分区** |
| VCE（方差分量） | `ps_vce` | ◐ | 秩 1 闭式化简；`tests/test_pipeline.py` 端到端使用并得到 `sig2_est`；**仍无逐量对拍** |
| 空间解缠（整数场 + 相位场） | `ps_spatial_unwrapping` | ◐ | 已修 `GAP-U01`（孤立节点零列导致法方程奇异、删弧后并行数组不同步）；`test_unwrapping_allows_isolated_node` + `test_pipeline` 端到端；仍缺 `fvl_pretest` 统计检验、删弧/分组迭代（当前校正循环在强噪声场景可能达上限） |
| 大气时间滤波 | `ps_filtering` | ✔ | `steps/filtering.py` + `pipeline.run_filtering`；写 `_psc_atmo_filt_sel*`，`tests/test_pipeline.py` 断言输出 |
| Kriging 大气估计 | `ps_atmo_kriging`/`ps_kriging`/`ps_fit_vario` | ◐ | `steps/kriging.py`：经验变差函数 + 拟合 + 普通 Kriging；未做 universal/局部块优化与端到端 |
| **PS 加密（首要需求）** | `ps_densification` | ◐ | `steps/densification.py` + `pipeline.run_densification`：KDTree 连接、批量 ILS、2of3/3of3，已写 `_psp_results_sel*`；端到端 65.6% 接受率、真 PS 富集（`test_densification_enriches_true_ps`）；仍无幅度分箱迭代、`mode_fix` |
| 形变模型（bowl） | `ps_defo_model`/`ps_defo_bowl` | ✘ | 参数表 `defo_method=[]`，默认不启用 |
| **噪声滤波（时间维平滑）** | `ps_noise_filter` | ✔ | `steps/noise_filter.py` + `pipeline.run_noise_filter`；`test_noise_filter_and_stc` 断言平滑降低时间粗糙度 |
| 时空一致性 | `ps_spatio_temporal_consistency` | ◐ | `steps/stc.py`：对 **defo（米）**做 STC（与 MATLAB 取数一致），邻域用 KD-tree 等价实现；仍未做网格缓冲分块 |
| 结果输出（`_defo`/`_trend`/`resfile`） | `ps_output` | ◐ | `pipeline.run_output` 写 `_psc_results_sel*`、`_ps_results_atmo_sel*`、`_psp_results_sel*`、点栈/CSV/resfile/summary；`master_res` 提供时写 `_psc_latlon_sel*` 与 CSV 经纬度；`_defo`/`_trend` 整幅栅格产物未实现 |
| 地理编码 | `geocoding_v0.8`/`metadataNew.m`/`intrp_orbit.m`/`lph2xyz.m` | ✔ | `io/doris.py` 读 `.res`+精密轨道，`core/orbit.py` 插值状态，`steps/geocoding.py` 实现 `lph2xyz/xyz2lp`；`test_doris_geocoding.py`（7 项）+ `run_output` 写 `_psc_latlon_sel*`。真实项目 `.res` 对拍待做 |
| 状态落盘（替代 `_project.mat`） | `save/load/pack` | ✘ | 每次 `depsipy run` 从头重算；只有 `_summary.json`，无跨进程续算 |

## 3. 参数语义待核（材料 vs 源码）

| 参数 | 材料描述 | 源码行为 | 处理 |
| --- | --- | --- | --- |
| `psc_selection_gridsize` | "网格大小"（表格给的是 m 单位） | `Nlines_g=round(grid/az_spacing)`，即**米** | 一致；预检按像素反算核对 |
| `defo_range=300` | 形变与大气相关的空间尺度 | `ps_vce`/`ps_kriging` 用作相关长度（m） | 待真实数据确认量级 |
| `ts_noise_filter_length=10/12` | "窗口长度"，源码注释为 `[year]` | `ps_noise_filter` 里按年→样本数换算 | 预检不自动换算，只提示 |
| `std_param(1)` | 表格写 "topo" | 源码单位是**米**，方差直接进 `Qy2` | 单位一致；数值是否合理需真实数据 |
| `defo_method=[]` | "未使用模型" | 第 6 组整段被跳过 | 预检报 `group6_skipped`（warning） |
| `max_mem_buffer=50e7` | "运行内存大小" | MATLAB 里是**元素数**上限，用于算缓冲 | Python 侧未使用该语义，待定 |
| `processing_groups=[1:9]` | "共 9 步" | 参考源码只有 1..8 | 预检 `step_mapping_pending`；基线启动器直接判定 blocked |

## 4. 端到端链路（2026-10-29，M2 MVP）

`depsipy run` 已把处理组 2..8 串成单进程链路（`src/depsipy/pipeline.py`）：

| 组 | 阶段 | 证据 |
| --- | --- | --- |
| 2 | selection + `_mrm/_amp_disp/_side_lobe_mask` + PSC/PSP 载荷 | `tests/test_pipeline.py::test_run_pipeline_writes_the_depsi_products` |
| 3 | VCE → 组网 → ILS → 解缠 → 每 PS 最终模型 | `test_step3_recovers_the_true_phase`（中位误差 <1 rad，0 次整数修正） |
| 3a | 去趋势（逐 IFG 平面），含"网络退化则回退"保护 | `test_detrend_keeps_the_phase_consistent`；`summary["detrend_accepted"]` |
| 4 | 二次组网 + 大气/非线性分离 | `test_run_pipeline...`（`_psc_atmo_filt_sel1.raw`） |
| 5/6 | 稠密化 | `test_densification_enriches_true_ps`（接受率、真 PS 富集、varfac） |
| 7/8 | STC、噪声滤波、结果落盘（含可选雷达地理编码） | `test_noise_filter_and_stc`、`test_output_geocodes_when_master_res_is_given` + 列布局断言 |

已知限制：解缠的整数一致性仍用启发式逐弧修正（强噪声下可能达到迭代上限，Matlab 侧同样
以 `Narcs*Nifgs` 为上限）；去趋势的拒绝循环加了"至少保留一半点"的保护（见 `docs/parity.md`
P-11），避免退化平面；加密是单趟（无幅度分箱与 `mode_fix`）；大气 Kriging 未接 PSP。

## 5. 本轮新增的可验证项

| 工具 | 作用 | 证据 |
| --- | --- | --- |
| `scripts/octave_parity.py` | 同一工程上 Octave（MATLAB 原码）vs depsipy 的 PSC/PSP/相位对拍 | `tests/test_matlab_parity.py`；本机 60×60×15：139/139 PSC、1456/1456 PSP 完全相同 |
| `scripts/readiness.py environment` | 记录解释器/依赖/线程/MATLAB/Octave | `tests/test_readiness.py` |
| `scripts/readiness.py preflight` | 只读校验参数、清单、栅格尺寸、样点 NaN | 15 项测试，含缺文件/尺寸/重复主影像/第 9 步 |
| `scripts/readiness.py compare` | 按坐标/日期身份对齐后比较，声明单位与参考点 | 覆盖重排、集合差、NaN/Inf、全 NaN |
| `scripts/readiness.py from-mat` | MATLAB `-v7` MAT → NPZ（仅显式字段） | `tests/test_interchange.py` |
| `scripts/readiness.py matlab` | 生成隔离的原生 MATLAB 基线计划（默认不执行） | 断言不加载 shim、不改原参数、输出重定向 |
| `scripts/matlab/capture_baseline.m` | 原生 MATLAB 采集：版本、Toolbox、函数解析、调用计数、产物 | 待真机执行（未验证） |
| `scripts/matlab/export_comparison.m` | 显式导出对照数组 | 待真机执行（未验证） |

## 6. 预检在本项目自身合成器上发现的缺陷（已修）

`simulate.generate_project` 曾把主影像**同时**写成一条干涉对和一条 `dummy` 行，导致
`orbit+date` 重复、`Nifgs` 与清单行数不自洽。预检的 `duplicate_acquisition` 直接命中。
现已改为「Nslc-1 条辅影像 + 1 条主影像」，并在 `truth.npz` 里显式存 `slave_index`。
合成器生成的项目现可通过预检（仅剩 warning/unverified）。

## 7. 待生产入口确认的问题

1. `step6`/`step9` 是否对应本仓库的 `group5/6` 与 `group8`，还是生产侧另有一层封装。
2. 生产环境是否在 `depsi.m` 之外还有前置（如 `create_psi_input_file.sh`）与后置脚本。
3. `dens_check='3of3'`、`Nest=3` 的实际判定与终止条件（参考源码 1282 行，含多分支）。
4. 结果产物（`_defo.raw`/`_trend.raw`/`resfile`）在当前生产中被哪些下游工具读取。
