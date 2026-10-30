# 架构

状态：**现行（2026-10-29）**。`io/`、`core/`、`steps/`、`simulate.py`、`pipeline.py`
与 `cli.py` 均已实现；§1 的 `runs/<id>/state.npz` 与 `manifest.json`（跨进程状态落盘）尚未实现，
其余数据流已按 `depsipy run` 串联并可执行。

## 1. 数据流

```text
param_file.txt ─┐
input_file.txt ─┴─► pipeline.py ─► [组1 初始化] ─► [组2 选择] ─► [组3 去趋势+组网]
                                        │              │              │
                                        ▼              ▼              ▼
                                  runs/<id>/       _mrm.raw      _dpsc_sel*.raw
                                  state.npz        _amp_disp.raw  _dpsc_results_sel*.raw
                                  manifest.json    _psc_2orig_sel*.raw
                                        │
                                        ▼
                     [组4 解缠+大气分离] ─► [组5 加密] ─► [组6 形变建模] ─► [组7 输出]
                              │                  │               │            │
                              ▼                  ▼               ▼            ▼
                       _psc_5defo_sel*.raw  _ps_4atmo.raw   _psc_defo_sel*  _defo.raw
                       _atmo_sel*.raw                                       _trend.raw
                                                                            _resfile.txt
                                                                            summary.json
```

## 2. 模块职责

| 模块 | 职责 | 对应 MATLAB |
| --- | --- | --- |
| `io/raw.py` | raw 栅格读写（memmap、1-based 裁剪、复数交错、稀疏像素） | `freadbk.m`/`fwritebk.m`/`freadcpx.m` |
| `io/params.py` | 参数文件解析与类型化访问 | `readinput.m`/`ps_readinput_parameters.m` |
| `core/models.py` | 27 个功能模型、断点、参数注解 | `ps_model_definitions.m` |
| `core/_reference.py` | LAMBDA 的纯 numpy 规范实现 | `ldldecom.m`/`decorrel_freek.m`/`ps_chistart.m`/`bk_lsearch_count.m` |
| `core/_kernels.py` | 上述的 numba 版本 + 批量入口 | 同上 |
| `core/lambda_ils.py` | 每模型的 LAMBDA 状态、逐弧/批量解析 | `ps_ils_and_bootstrap_estimation.m` |
| `core/estimation.py` | 几何（`Qahat`/`rhs_fixed`/OMT/相干性）、多模型选择 | 同上 + `ps_vce.m` |
| `steps/selection.py` | 幅度流式统计、PSC/PSP、旁瓣掩膜、网格选取 | `ps_selection.m` |
| `steps/network.py` | Delaunay/spider 弧网、独立弧、方差分量 | `ps_form_network.m`/`ps_vce.m` |
| `steps/unwrapping.py` | 整数场最小二乘 + 修正，相位最小二乘 | `ps_spatial_unwrapping.m` |
| `io/doris.py` | Doris/RapidSAR 产物读取（MFF 头、`.r00/.x00/.j00`、`.res` 元数据/轨道） | `metadataNew.m`（Doris `.res` 约定） |
| `core/orbit.py` | 精密轨道多项式拟合 + 样条状态插值 | `intrp_orbit.m` |
| `steps/geocoding.py` | WGS84↔ECEF、`lph2xyz`/`xyz2t`/`xyz2lp`、雷达→大地 | `lph2xyz.m`/`xyz2t.m`/`xyz2lp.m`/`xyz2ell.m` |
| `steps/densification.py` | KDTree 连接 + 批量 ILS + 2of3/3of3 | `ps_densification.m` |
| `steps/filtering.py` `kriging.py` | 大气/非线性分离、变差函数与 Kriging | `ps_filtering.m`/`ps_atmo_kriging.m` |
| `steps/noise_filter.py` `stc.py` | 时间维平滑、时空一致性 | `ps_noise_filter.m`/`ps_spatio_temporal_consistency.m` |
| `pipeline.py` `cli.py` | 处理组 2..8 编排与 `depsipy run` | `depsi.m` |
| `simulate.py` | 合成工程 + 真值 | `depsi_sim.m`（重构版） |

## 3. 状态格式（**未实现**，M2 只剩这一项）

当前每次 `depsipy run` 从头重算，只写 `_summary.json`。规划的跨进程续算格式如下（尚未落地）：

`runs/<project_id>/manifest.json`

```json
{
  "schema": 1,
  "project_id": "sim",
  "sensor": "spacety",
  "Nlines": 120, "Npixels": 120,
  "Nslc": 12, "Nifgs": 11,
  "crop": [1, 120, 1, 120],
  "groups_done": [1, 2, 3],
  "timings": {"1": 0.4, "2": 12.1, "3": 30.5},
  "params_snapshot": {"ps_threshold": 0.3, "...": "..."}
}
```

`runs/<project_id>/state.npz`：数组（`sig2_est`、`psc_*`、`arcs`、`acheck_*`…），键名与 MATLAB 变量同名，便于对读。

## 4. 与 MATLAB 版本的文件兼容

写出同名 `.raw`，布局逐字节一致（行主序 + 复像素交错 + double/float32）。因此：

- 现有 DePSI 工程可以直接被 `depsipy` 处理（无需迁移数据）；
- `depsipy` 的产物也能被 `ps_output.m` 读取（M6 核对）。

差异记录见 [parity.md](parity.md)。
