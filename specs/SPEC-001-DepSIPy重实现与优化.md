# SPEC-001 · DePSI 的 Python 重实现与性能优化

状态：**草案 / 待 Review**（会话式迭代第 1 轮）
作者：agent
日期：2026-10-29

---

## 1. 背景与事实（已核实）

工作区 `workspace/` 下有三个来源：

| 目录 | 性质 | 事实 |
| --- | --- | --- |
| `ps_toolbox_v1.7.8.8_hhf/` | MATLAB · 主线 | DePSI v1.7.8.8（TU Delft PS-InSAR 时序列处理），`main/` 48 个 `.m`、`simulation/` 20 个、`lambda/` 19 个（LAMBDA 整数最小二乘）、`quadtree/`、`geocoding_v0.8/`。`main/depsi.m` 1403 行，分 1–8 处理组。本工作区版本被二次修改过：`ps_ils_and_bootstrap_estimation.m` 里加了 `parfor` + CPU 自动调度，`ps_selection.m` 里加了 GexYuu 旁瓣抑制，`depsi.m` 里加了 `g_log_file` 计时。 |
| `ps_toolbox/ps_toolbox/` | MATLAB · 旧线 | v1.7.7.5，与主线同构，仅 `depsi.m` 1196 行（旧版）。可视为参考/回退点。 |
| `Doris干涉加时序depsi/doris干涉/` | Python | DORIS 干涉处理脚本（RapidSAR，141 个 `.py`）：配准、干涉图、去平地、滤波、解缠、地理编码。**不含时序（PS）反演**。 |
| `Doris干涉加时序depsi/depsi时序/` | 胶水 | `depsi_test.py` 只是 `os.system("matlab -nodisplay -r run('depsi_test.m')")`。即：**时序部分目前完全依赖 MATLAB**。 |
| `InSAR_affair-shallow-2026-10-19/` | 独立项目 | DANI-NET 论文复现（PyTorch/昇腾），与本任务无关，仅共享 InSAR 领域。 |

结论：用户所说的「原有 matlab 程序速度慢，python 版本功能不全」，准确含义是
**干涉链已有 Python（RapidSAR），时序列（DePSI）只有 MATLAB，且慢**。

### 1.1 已确认的瓶颈（来自代码事实）

1. `ps_selection.m`：对每个网格缓冲区把 **整块 slc/ifg/h2ph 全部读入**（`freadbk` 逐行 `fseek+fread`），而真正用到的只有候选像素。I/O 量 ≈ 数据总量 × 常数，与候选数无关。
2. `ps_ils_and_bootstrap_estimation.m`：每弧调用 `inv(C'*invQy*C)`、`inv(Qy)`、`decorrel`，再逐弧跑 MATLAB 标量循环的 LAMBDA 搜索；`parfor` 只能缓解不能改变单弧常数因子过大。
3. `ps_vce.m`：`N(k,j)=trace(QP·Qy1_k·QP·Qy1_j)` 直接按定义算 → `Nsig²` 次 `Nifgs×Nifgs` 矩阵乘，复杂度 O(Nsig²·Nifgs³)。
4. `ps_spatial_unwrapping.m`：500+ 行逐弧/逐 ifg 循环 + 每轮重算 `chol`、`inv(Qyy)`。
5. 全局 `pack` / `save/load *_project.mat` 在每组之间整工作区落盘。
6. 数值风险：`inv(C'*invQy*C)` 在高度-相位因子近似常数时病态（cond ≈ 1e10 量级），会产出 1e15 量级参数；实测本仓库代码路径确实会触发。

### 1.2 已确认的语义细节（重实现必须对齐）

- raw 格式：**行主序 + 复像素交错**，等于 numpy `complex64` 的 C 序，无头；宽度由文件大小反推。
- 参数文件：`textread('%s')` + `str2num`（即 `eval`），支持 `[1:9]`、`[30,100]`、`[]`、带引号字符串。
- `input_file` 列序：`orbit date slc ifg h2ph out Bperp Btemp Bdop l0 lN p0 pN`，主影像行用 `dummy`。
- 模型库 27 个（`ps_model_definitions.m`），参数 1..10：topo / master / subpixel / linear1-3 / quadratic / cubic / periodic1-2。
- 估计器：`Qahat = inv(C'invQy C)[1:Nifgs]`，LAMBDA 去相关 + bootstrap + 分支定界，`varfac = OMT/b`（**用不带先验的 `rhs_fixed` 拟合**，与搜索用的正则化度量不同——原文如此）。
- 解缠：先对**整数模糊度场**做加权最小二乘 + 逐弧修正，再对**缠绕相位**做同一法方程求解，`unw = 2π·acheck + φ_smooth`。

---

## 2. 目标 / 非目标

**目标**

1. 交付 `depsipy`：可直接读现有 `param_file.txt` / `input_file.txt` / `*.raw` 的 Python 时序列处理器，功能覆盖到 DePSI 的 1–8 组（选择、网络、估计、解缠、分离、加密、形变建模、输出、噪声滤波）。
2. 性能：同一数据规模下，端到端显著快于 MATLAB（目标 ≥ 10×，热点步骤 ≥ 30×），内存可控（流式 I/O，不整栈入内存）。
3. 可验证：自带合成数据生成器，用**已知真值**验收（形变速率/地形/时序 RMS），并有单元测试锁定每个数值核。
4. 与原实现可对照：文件布局、命名、resfile 尽量兼容；差异逐条写入 `docs/parity.md`。

**非目标（本轮）**

- 不改 RapidSAR 干涉链（另一条线，只在接口处对齐）。
- 不重写 `simulation/depsi_sim.m` 的全部统计评估（只提供足以验收的合成器）。
- 不做 GUI / 云端 / 数据库。
- 不追求与 MATLAB 逐位相同；追求**同一数学模型下的等价结果**，差异需有解释。

---

## 3. 约束

- 运行环境：`python3.13`，`numpy 2.5` / `scipy 1.18` / `numba 0.68`（可选），16 核 / 14 GB / 无 GPU；venv 在 `workspace/depsipy/.venv`。
- 无 MATLAB、无 Octave，无法在线对照运行；因此**用合成器 + 数学不变量**代替对拍（见 §5）。
- 许可证：原工具箱声明 "scientific use only, commercial prohibited"。`depsipy` 是独立重实现，但算法与文件格式沿用；商用需用户另行确认授权路径（**待用户决策**）。
- 会话式迭代：每轮汇报「完成了什么 / 证据 / 缺什么 / 下一步需谁确认」，由用户裁决优先级。

---

## 4. 方案与取舍

| 决策点 | 选择 | 备选 | 理由 |
| --- | --- | --- | --- |
| 语言/栈 | numpy + scipy(+numba) | C++/Rust 扩展、Julia | 生态与可维护性；numba 足以消除标量循环 |
| 数值核 | 逐弧向量化 + numba 并行 | 纯 numpy、Cython | 单弧 LAMBDA 是串行分支定界，只能编译化 |
| I/O | `memmap` + 行块流式 + 稀疏像素读取 | 整栈读入 | 与 MATLAB 的差距主要来自 I/O |
| 状态 | `npz` + `json`（替代 `_project.mat`） | HDF5 | 零依赖；大数组仍可 mmap |
| 解缠 | 同一法方程、稀疏因子化一次解全部 ifg | 逐弧循环 | 数学等价、量级提速 |
| VCE | 用 `Qy1_k` 的秩 1 结构把 trace 化简为闭式 | 按定义算 trace | 等价且把 O(Nsig²·Nifgs³) 降为 O(Nsig²) |
| 验收 | 合成数据 + 真值 RMS + 单元测试 | 与 MATLAB 逐位对拍 | 环境无 MATLAB；且对拍本身不验证物理正确性 |

---

## 5. 任务拆解与当前状态

| # | 任务 | 状态 | 证据 |
| --- | --- | --- | --- |
| T1 | raw I/O（`cpxfloat32`/`cpxint16`/裁剪/稀疏像素） | **已完成·已验证** | `tests/test_io_raw.py` 8 项 |
| T2 | 参数文件解析（范围/多行/引号） | **已完成·已验证** | `tests/test_io_params.py` 9 项，含真实 `param_file.txt` |
| T3 | 27 个功能模型 + 断点模型 | **已完成·已验证** | `test_functional_model_columns`、`test_breakpoint_model_*` |
| T4 | LAMBDA（`ldldecom`/`decorrel`/`chistart`/`lsearch`）+ numba 核 | **已完成·已验证** | 与穷举枚举对拍；numba 与参考实现 12 组参数一致 |
| T5 | 几何/估计（`Qahat`、`rhs_fixed`、OMT、ensemble coherence） | **已完成·已验证** | 噪声自由时参数精确复原；ILS 在 `Qahat` 度量下对 ±1 邻域最优 |
| T6 | PSC/PSP 选择 + 旁瓣掩膜 + 网格选取 | **代码完成·待验证** | `steps/selection.py`，尚无端到端证据 |
| T7 | 弧网（Delaunay / spider）+ VCE | **代码完成·部分验证** | VCE 秩 1 化简待单测 |
| T8 | 空间解缠（整数场 LS + 相位 LS） | **代码完成·已验证核心** | 合成 PS 场：整数模糊度正确时误差 0；20/407 弧被破坏仍可完全恢复 |
| T9 | 合成数据生成器（真值） | **代码完成·待跑通** | `simulate.py` |
| T10 | 每组串起来的 pipeline + CLI + 状态落盘 | **进行中** | `pipeline.py`/`cli.py` 未写 |
| T11 | 分离大气（时间滤波 + Kriging） | **未开始** | — |
| T12 | PS 加密、形变模型（bowl）、噪声滤波、STC | **未开始** | — |
| T13 | 输出（`_defo/_trend/_resfile`、GeoTIFF/CSV）与可视化 | **未开始** | — |
| T14 | 基准报告（vs MATLAB 复杂度/实测）与 `docs/parity.md` | **未开始** | — |

---

## 6. 验收标准（本轮要达成的）

1. `depsipy run <param_file>` 在合成项目上跑完 1→7 组，产出形变/速率图与 resfile。
2. 合成数据（Nlines=Npixels=120、Nslc=12、PS 占比 ~8%）上：
   - PSC 找回率 ≥ 80%（真 PS 中位幅度离散 < 0.2 的部分）；
   - 主分量 PS 的**线性速率**与真值相关系数 ≥ 0.95，RMS 误差 ≤ 30% 真值标准差；
   - 参考点、PS 数、弧数写进 resfile。
3. 单元测试全绿（当前 48 项）。
4. 性能基准：`benchmarks/` 下给出各步骤耗时表；选择与估计两步相对 MATLAB 同算法的理论复杂度下降 ≥ 10×（用 I/O 字节数与标量运算数估算，并标注为估算）。

---

## 7. 风险

| 风险 | 影响 | 处置 |
| --- | --- | --- |
| 无 MATLAB 对照 | 语义偏差不可见 | 逐函数对读源码 + 记录 `docs/parity.md`；关键处写不变量测试 |
| `Qahat` 病态（h2ph 近似常数） | 参数爆量 | 已改为 solve/Cholesky + 条件数检查 + pinv 兜底并告警 |
| 用户真实数据巨大（GB 级 raw） | 内存/时间 | 流式 + memmap；`max_mem_buffer` 语义改为字节预算 |
| 断点模型/加密/形变模型等复杂分支 | 功能缺口 | 按 §5 优先级逐轮补齐，缺口显式登记 |
| 许可证 | 商用受限 | 待用户决策，不擅自改声明 |

---

## 8. 验证命令

```bash
cd workspace/depsipy
.venv/bin/python -m pytest tests -q                     # 单元测试
PYTHONPATH=src .venv/bin/python -m depsipy.cli simulate --out /tmp/sim1   # 生成合成项目
PYTHONPATH=src .venv/bin/python -m depsipy.cli run /tmp/sim1/param_file.txt
PYTHONPATH=src .venv/bin/python -m depsipy.cli validate /tmp/sim1         # 与真值比对
```

## 9. 待用户决策

1. 是否在 `workspace/depsipy/` 初始化 `AGENTS.md` + `rules/` + `specs/` + `docs/` 治理结构（本文件所在位置即为 `specs/`）。
2. 下一轮优先级：**(a)** 打通 1→7 组端到端（推荐，先拿到可用产物）；**(b)** 先把选择/估计做到极致性能并出基准；**(c)** 先补 11–12 组（大气分离/加密/形变模型）再打通。
3. 商用授权口径（是否允许沿用 DePSI 算法与文件格式对外分发）。
