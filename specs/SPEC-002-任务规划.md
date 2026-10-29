# SPEC-002 · 任务规划（里程碑）

状态：**建议 / 待用户确认**
日期：2026-10-29
前置：`ADR-001-技术栈.md`（栈选型）、`SPEC-001-DepSIPy重实现与优化.md`（背景与瓶颈）

工作方式：**会话式迭代**。每个里程碑结束都产出「完成 / 证据 / 缺口 / 待决」，由用户裁决是否进入下一里程碑。

---

## 里程碑总览

| M | 名称 | 交付物 | 状态 |
| --- | --- | --- | --- |
| M0 | 数值地基 | raw I/O、参数解析、27 模型、LAMBDA/ILS、估计几何、48 项测试 | **已完成** |
| M1 | 合成数据 + 选择 | 合成项目生成器、PSC/PSP 选择、mrm/amp_disp/旁瓣图、选择验收 | **已验收**（见 §M1 实测） |
| M2 | 端到端主链 | pipeline（1→7 组骨架）、CLI、状态落盘、速率/形变输出、真值验收 | 进行中 |
| M3 | 大气与形变分离 | 时间滤波 + Kriging（第 4 组），形变时序产物 | 未开始 |
| M4 | 性能工程与基准 | 流式/批量改造、内存剖析、`benchmarks/report.md`、10× 目标验证 | 未开始 |
| M5 | 加密与时序精化 | PS 加密、形变模型（bowl）、噪声滤波、时空一致性（第 5/6/8 组） | 未开始 |
| M6 | 输出与互通 | GeoTIFF/CSV/JSON、与 RapidSAR 的串联脚本、MATLAB 兼容性核对 | 未开始 |
| M7 | 文档与交接 | `docs/parity.md` 完整、用户手册、示例工程、复现脚本 | 未开始 |

依赖：M1 → M2 → M3 → M5；M4 可在 M2 后任意时点插入；M6 依赖 M3；M7 与各 M 并行推进。

---

## M0（已完成）

- 证据：`.venv/bin/python -m pytest tests -q` → 48 passed。
- 关键数字：LAMBDA numba 相对纯 Python **37×**（97k 弧/s，n=10）。
- 已知缺口：`steps/` 下的选择、网络、解缠尚未有端到端证据。

## M1 · 合成数据 + 选择

**交付**
1. `simulate.generate_project()` 能生成 DePSI 格式工程 + `truth.npz`（已写，待跑通）。
2. `steps/selection.select_psc_psp()` 输出 PSC/PSP 列表 + `_mrm.raw`/`_amp_disp.raw`/`_side_lobe_mask.raw`（已写）。
3. 选择验收脚本：与 `truth.ps_mask` 比对。

**验收**
- 合成工程（120×120、12 景、PS ~10%、DS 幅度调制 σ=0.45、PS σ=0.12）：
  - **精度**（选出的 PSC 中真 PS 占比）≥ 95%；
  - **网格单元覆盖率**（含可达真 PS 的单元产出 PSC 的比例）≥ 95%；
  - `_amp_disp`：真 PS 中位 < 0.2，真 DS 中位 > 0.3；
  - 旁瓣掩膜误杀真 PS < 5%。
  > 实测（seed=0，120×120，12 景）：精度 1.00、覆盖率 1.00、PS 中位 0.117 / DS 0.426、
  > 旁瓣误杀 0.00、耗时 0.03 s。详见 `docs/selection-calibration.md`。
- 选择耗时随 raster 面积线性、随候选数近似不变（给出 I/O 字节数对比 MATLAB 的估算）。

**风险**：旁瓣掩膜是原作者的启发式，参数（`peak_tolerance=0.7`、95 分位）在合成数据上需重新校准；校准过程写进文档。

## M2 · 端到端主链（1→7 组骨架）

**交付**
1. `pipeline.py`：读参数 → 几何常量 → 组 1..7 顺序执行，状态存 `runs/<id>/state.npz` + `manifest.json`。
2. `cli.py`：`depsipy run|simulate|validate`（argparse，无额外依赖）。
3. 输出：`_defo.raw`（各景形变，mm）、`_trend.raw`（速率，mm/yr）、`_psc_results_sel*`、`_resfile.txt`、`summary.json`。
4. 每 PS 的最终参数估计（`final_model`，默认 2）与 `ens_coh`、`sig2hat`。

**验收**（合成数据，同 M1 配置）
- 主分量 PS 的速率与真值：`Pearson r ≥ 0.95`，`RMS ≤ 30% × std(truth)`；
- 地形：相关系数 ≥ 0.9（有断层/断点时允许分段）；
- `ens_coh` 中位数 ≥ 0.8；
- resfile 记录 PS 数、弧数、参考点、各组耗时。

**风险**：`Qahat` 在 `h2ph` 近常数时病态 → 已加条件数告警与 pinv 兜底；若真数据出现大面积参数爆量，需要重新审视 `std_param` 默认值（原默认 `[30,...,0.005]` 对空间基线稀疏的数据并不合理）。

## M3 · 大气与形变分离（第 4 组）

**交付**
1. `steps/filtering.py`：时间域高斯滤波分离大气与非线性形变（`ts_atmo_filter_length`）。
2. `steps/kriging.py`：变差函数拟合（`ps_variogram`/`ps_fit_vario`）+ 普通 Kriging 插值大气相位，输出 `_atmo_sel*`、`_defo_sel*`。
3. 重新组网 + 二次解缠（`results_id='atmo'`，对应 `depsi.m` 组 4 的第二段）。

**验收**
- 合成数据（注入已知 APS：空间相关长度 ~150 m，幅度 0.6 rad）：分离后 APS 残差标准差较分离前下降 ≥ 50%；
- 非线性形变（季节性 1 cm）恢复误差 RMS ≤ 40% 幅度；
- Kriging 预测与真值 APS 的相关系数 ≥ 0.7。

**风险**：Kriging 是 O(N³)（每 IFG 一次）→ 需要分块/邻域截断；若 PS 数 > 5×10⁴ 必须做局部 Kriging（写进实现）。

## M4 · 性能工程与基准

**交付**
1. `benchmarks/run_bench.py`：规模参数化（small/medium/large），输出耗时、峰值内存、I/O 字节数。
2. 优化项：选择步骤的块大小自适应；估计的批处理与缓存复用；解缠的因子化复用；避免重复读 raw。
3. `benchmarks/out/report.md`：与 MATLAB 的**估算**对比（复杂度 + I/O 字节），并标注「无 MATLAB 环境，未实测」。

**验收**
- 端到端相对 M2 首个可用版本 ≥ 3×；相对 MATLAB 的估算 ≥ 10×（需在报告里给推导）。
- 峰值内存 ≤ 2 GB（120×120×12 的合成项目），且随 raster 面积线性而非平方增长。

## M5 · 加密与时序精化（第 5/6/8 组）

**交付**：PS 加密（幅度离散分箱 + 迭代）、形变模型（bowl: `ps_defo_bowl`）、噪声滤波（时间域高斯）、时空一致性指标（STC）。
**验收**：加密后 PS 数 ≥ 2× 且新增 PS 的 `ens_coh` 中位数 ≥ 0.7；STC 通过率与 MATLAB 同量级（无 MATLAB 时用分布对比并记录）。
**风险**：加密模块 1282 行、分支多，是剩余工作里最重的一块；建议单独一轮迭代。

## M6 · 输出与互通

**交付**：GeoTIFF（可选 `rasterio`）、CSV 点位表、JSON 摘要；`depsipy run-rapidsar`（或 shell 脚本）把 RapidSAR 的干涉产物直接喂进时序列；MATLAB 兼容性核对清单（同名 `.raw` 能被 `ps_output.m` 读）。
**验收**：`ps_output.m` 的读取路径能在我们的产物上跑通（无 MATLAB 时，用逐字节格式核对替代并说明）。

## M7 · 文档与交接

**交付**：`docs/parity.md`（每条差异：MATLAB 行为 / depsipy 行为 / 原因 / 影响 / 开关）、用户手册、示例工程（`examples/sim_small/`）、一键复现脚本。

---

## 本轮（M1）的具体动作清单

1. 跑通 `simulate`，生成 `/tmp/dp_sim`（含 truth）。
2. 写 `tests/test_simulation.py`：验证生成文件可被 `io/raw` 读回、truth 自洽（形变/地形量级）。
3. 写 `tests/test_selection.py`：在选择结果上算 PSC 找回率/假阳性，锁定 M1 验收门槛。
4. 若门槛不达标：调 `psc_threshold`、旁瓣参数、`ps_fraction` 的合成设定，把校准过程与结论写入 `docs/selection-calibration.md`。
5. 汇报：完成/证据/缺口/下一步。

---

## 待用户确认

1. 是否认可 `ADR-001` 的栈选型（Python + numpy/scipy + numba，raw 兼容，npz 状态）。
2. 里程碑顺序是否认可（M1 → M2 → M3 → M5，M4 可提前插入）。
3. M3 的验收里「APS 残差下降 ≥50%」是否是你关心的口径，还是更关心**时序形变曲线**的保真度。
