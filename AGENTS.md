# AGENTS.md · depsipy

本文件是本项目的**长期记忆**：已核实的事实、可复现命令、边界与当前决策。
未核实的内容一律标注「待确认」，不得当作规则使用。

- 项目根：`/home/zoneip/工程/Zoneip/depsi`（2026-10-29 已从 `workspace/depsipy` 展开到根目录）
- 参考工作区：`workspace/`（含 MATLAB 原实现、Doris Python 链；不纳入仓库）
- 工作方式：**会话式迭代 + 用户主导**（用户给方向与顺序；每轮汇报「完成/证据/缺口/待决」）
  - 用户 2026-10-29 指示：这是长任务，存在相当多不确定决策，**按用户思路推进，不做过多无用功**；
    数据与服务器随后提供（用户称可 100% 跑大规模任务）。因此：**不预先实现未经确认的功能**，
    不跑大规模合成基准，不为了"完整"而扩写文档；先等输入再动手。

## 1. 这是什么

`depsipy` 是 TU Delft **DePSI v1.7.8.8**（Persistent Scatterer InSAR 时序列处理）
的 Python 重实现。目标是：读同一套 `param_file.txt` / `input_file.txt` / `*.raw`，
跑同一条处理链，但显著更快、内存可控、可验证，并补齐原 Python 侧缺失的功能。

上游事实：
- MATLAB 主线：`workspace/ps_toolbox_v1.7.8.8_hhf/`（被二次修改，见 SPEC-001 §1）
- MATLAB 旧线：`workspace/ps_toolbox/ps_toolbox/`（v1.7.7.5，参考/回退）
- 干涉链（Python，RapidSAR/DORIS）：`workspace/Doris干涉加时序depsi/doris干涉/`
- 时序侧原先靠 `depsi_test.py` 调 `matlab -nodisplay`；**本项目要替换的就是这一步**

## 2. 目录职责

```text
src/depsipy/
  io/        raw 栅格（memmap、1-based 裁剪、复数交错）、参数文件、工程状态
  core/      功能模型、LAMBDA/ILS、估计几何（Qahat/rhs_fixed/OMT）、编译核
  steps/     处理组：selection / network / unwrapping / filtering / ...
  simulate.py  合成项目生成器（写 DePSI 格式 + truth.npz）
  cli.py       `depsipy run|simulate|validate`
tests/       单元测试（pytest）
reference/   MATLAB 源码只读副本（本地对读用，不纳入公开仓库）
specs/       方案、ADR、验收标准
rules/       现行架构/编码/质量门规则
docs/        架构说明、对拍记录、用户文档
achieve/     归档（只归档，不删除）
benchmarks/  耗时与规模基准
```

依赖方向：`steps → core → io`，`core` 不依赖 `steps`；`io` 不依赖数值模型。
`reference/` 与 `simulate.py` 不得被 `src` 内部 import（仅测试/CLI 使用）。

## 3. 可复现命令

```bash
cd depsi
python3 -m venv .venv && .venv/bin/pip install -e '.[fast,plot]'   # 首次
.venv/bin/python -m pytest tests -q                                 # 质量门 1
PYTHONPATH=src .venv/bin/python -m depsipy.cli --help
PYTHONPATH=src .venv/bin/python benchmarks/run_bench.py             # 基准
```

环境（2026-10-29 实测）：Python 3.13.5、numpy 2.5.4、scipy 1.18.1、numba 0.68.1、
16 vCPU、14 GiB RAM、**无 GPU**、**无 MATLAB**；已装 **Octave 9.4.0** 仅用于对拍（§4.8）。

算力规划（用户 2026-10-29 告知）：后续会提供 **64 × 2 GB 显存的大机器**，
算力不再是瓶颈；届时优先把向量化/批量核迁到该机器跑大规模验证。具体型号、
CPU/内存/磁盘、驱动与 CUDA/加速器栈**待确认**；在拿到机器前仍按本机规格设计。

## 4. 已核实的关键事实（改动前必读）

1. raw 布局 = 行主序 + 复像素交错 + 无头；宽度由文件大小反推（`io/raw.py`）。
2. 参数文件是 `textread('%s')` + `str2num`（本质 `eval`）：支持 `[1:9]`、`[30,100]`、`[]`、引号串。
3. `input_file` 列序：`orbit date slc ifg h2ph out Bperp Btemp Bdop l0 lN p0 pN`，主影像行 `dummy`。
4. 功能模型 27 个（`core/models.py`），参数编号 1..10 与 MATLAB 一致。
5. `varfac` 用**不带先验**的 `rhs_fixed` 拟合；而 ILS 搜索用**带先验**的 `Qahat` 度量。
   两者不是同一统计量，测试不可互相对拍（见 `docs/parity.md`）。
6. 解缠 = 先整数场最小二乘 + 逐弧修正，再相位最小二乘，`unw = 2π·acheck + φ_smooth`。
7. MATLAB 的 `max_backtrack = n³` 会让 LAMBDA 提前放弃并回退到 bootstrap（非最优）；
   `depsipy` 默认放宽该上限，`compat=True` 可还原。
8. 对拍通道：本机已装 **Octave 9.4**（含 image/statistics/parallel），配合
   `benchmarks/octave/shim/` 的 7 个 shim 可跑 MATLAB 原码。注意两点：
   (a) `ps_set_globals` 会把同名**函数参数**遮蔽成空全局，驱动脚本必须先赋值；
   (b) Octave 的 `parfor` 串行，故 MATLAB 侧耗时是"未并行"的下界。

## 5. 分支规划

- `master`：当前 Python 实现（主开发线）。
- `rust`：预留分支，未来 Rust 版本在此线开发；创建于 2026-10-29，与 master 同步起点。
- 语言无关的规格（`specs/`）、规则（`rules/`）、对拍记录（`docs/parity.md`）与基准脚本
  保留在 master；Rust 分支按需同步，不单独复制一份。

## 6. 当前状态（每轮更新）

- 已完成·已验证：raw I/O、参数解析、27 模型、LAMBDA+ILS、估计几何、合成器、
  PSC/PSP 选择（精度 1.00 / 单元覆盖率 1.00，见 docs/selection-calibration.md）；56 项测试
- 代码完成·部分验证：弧网+VCE（VCE 待单测）、空间解缠（合成场已验证核心）
- 进行中：pipeline/CLI/状态落盘（M2）
- **性能对拍已完成（v0）**：estimate 987×（3000 弧×19 ifg，Octave 10.07s → 0.0102s，
  整数解 3000/3000 相同）、selection 10.8×（400×400×20，PSC 逐个相同）；详见
  `docs/performance.md` 与 `benchmarks/results/`。对拍还发现并修复 2 个真 bug（见 §4.8）。
- 未开始：大气分离/Kriging、加密、形变模型、噪声滤波、STC、输出与可视化、基准
- 授权：沿用「仅限科研，禁止商用」；2026-10-29 用户决策建立公开仓库（代码公开，许可不变）。

## 7. 纪律

- 改数值核前先跑 `pytest`；新增分支要有对应测试或不变量。
- 不与 MATLAB 逐位对齐为目标；对齐的是**数学模型**，差异写入 `docs/parity.md`。
- 不提交 `.venv/`、`*.raw`、`truth.npz`、`runs/`。
- 不擅自提交/推送/删除归档；需要时先问。

## 8. 待用户提供（阻塞项）

| # | 需要什么 | 用来决定什么 | 没它时不做的事 |
| --- | --- | --- | --- |
| 1 | 一份真实工程：`param_file.txt` + `input_file.txt` + 若干 `*.raw`（可只读） | 真实幅度/相位分布、`std_param` 是否合理、预处理链（RapidSAR）产物是否直接可读 | 调选择阈值、定性能常数、判断老版精度问题是否复现 |
| 2 | ~~服务器画像~~ 已部分告知：未来提供 64 × 2 GB 显存大机器，算力不受限；具体型号/CPU/内存/磁盘/加速栈待确认 | 并行策略、是否上 GPU 核 | 在拿到机器前不做 GPU 专用实现 |
| 3 | 优先级：先"跑通出产物"还是先"补功能（加密/形变模型/大气分离）" | M2 与 M5 的顺序 | 继续 M2 之后的实现 |

已知可立即交付的东西（不需要新信息）：`depsipy` 的 raw I/O、参数解析、模型、LAMBDA/ILS、
选择、弧网、解缠；CLI 与 pipeline 尚未落地（`src/depsipy/cli.py`、`pipeline.py` 不存在）。
