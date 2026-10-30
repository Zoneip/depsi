# AGENTS.md · depsipy

本文件是本项目的**长期记忆**：已核实的事实、可复现命令、边界与当前决策。
未核实的内容一律标注「待确认」，不得当作规则使用。
你的预算将是无限的，需要任何算力、设备、参考都可以停下来找我！！！

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
- 官方 Python 版：`reference/python/` = TUDelftGeodesy/DePSI（Apache-2.0，beta，
  commit `51939cb`，2026-10-28）。架构与本项目不同：xarray/Dask + 模块化（network /
  arc_estimation / atmosphere_estimation / densification / point_quality），不含 MATLAB 的
  1–8 处理组一一对应关系；其 MATLAB 版在该仓库 `depsi-matlab` 分支（仓库 README 写的是
  `stable`，该分支现已不存在）。
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
scripts/     上云前工具：readiness.py + matlab/ 采集脚本
reference/   参考源码（仅本地，不纳入公开仓库）
  matlab/      TU Delft MATLAB v1.7.8.8 原码（对读用）
  python/      官方 Python 版 TUDelftGeodesy/DePSI（Apache-2.0，对读/架构参考）
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
.venv/bin/python -m pytest tests -q                                 # 质量门 1（124 passed；云端 122 + 2 skipped）
.venv/bin/python scripts/readiness.py --help                        # 上云前准备工具
.venv/bin/python scripts/readiness.py environment
.venv/bin/python scripts/readiness.py preflight <param_file.txt> [--map pathmap.json]
.venv/bin/python scripts/readiness.py matlab --entry <depsi.m> --params <param> --run-dir runs/base1
PYTHONPATH=src .venv/bin/python benchmarks/run_bench.py             # 基准
```

注意：`depsipy run` 已可用（`src/depsipy/pipeline.py` + `cli.py`），按
`processing_groups` 串起 selection → 组网/估计/解缠 → 去趋势 → 加密 → 大气滤波 →
噪声/STC → 输出；`pyproject.toml` 声明 `depsipy = depsipy.cli:main`。

环境（2026-10-29 实测）：Python 3.13.5、numpy 2.5.4、scipy 1.18.1、numba 0.68.1、
16 vCPU、14 GiB RAM、**本机当前无 GPU**、**本机当前无 MATLAB**；已装 **Octave 9.4.0** 仅用于对拍（§4.8）。

MATLAB 环境（用户 2026-10-29 告知）：**可 100% 提供 MATLAB 环境**；版本、Toolbox
（尤其 Parallel Computing Toolbox）、编译器与接入时间**待确认**。在 MATLAB 接入前，
不以“本机无 MATLAB”为由否决 MEX 等方案；接入后用真机实测再决策。

算力规划（用户 2026-10-29 告知）：后续会提供 **64 × 2 GB 显存的大机器**，
算力不再是瓶颈；届时优先把向量化/批量核迁到该机器跑大规模验证。具体型号、
CPU/内存/磁盘、驱动与 CUDA/加速器栈**待确认**；在拿到机器前仍按本机规格设计。

实测服务器（2026-10-29 临时接入；同日 19:17 再次接入复核）：aarch64（HiSilicon Kunpeng，
640 vCPU / 8 NUMA / 1.9 TiB RAM，load≈34 属**共享机**）、Ubuntu 22.04.5 容器（宿主内核
openEuler 5.10）、CANN 9.0.0 在 `/usr/local/Ascend`、overlay 磁盘 300 GB（剩 56 GB）、
`/workspace` 为 glusterfs 18 TB 共享盘、Python 3.10.12（系统级，**无 pip3**）、pytest 9.2.0、
gcc 11.4 + OpenMP。`nvidia-smi` 不存在；`npu-smi` 因缺 `libc_sec.so` 无法启动（待修）。
外网：pypi 200 可达，**github 不可达**（同步只能走 SSH 通道，不能 git clone）。
**该节点是 ARM64，官方 MATLAB Linux 版仅支持 Intel/AMD x86-64，因此 MATLAB 与 MEX 不能装在
该 ARM 节点**；其算力可用于 Python/NPU 或原生 ARM 代码（Octave/Rust）。
工程在 `/workspace/depsi`（无 `.git`，按文件同步）；其 `.venv` 是**伪 venv**
（`.venv/bin/python` → `/usr/bin/python3`，依赖实际来自系统 site-packages：
numpy 2.2.6 / scipy 1.15.3 / numba 0.68.1）。
2026-10-29 19:20 整树重同步后 pytest 收集 110 项、跑到 100%、**无 failure/error**
（进度行可见 1 skip + xfail 标记；确切计数行未捕获，因连接随即被平台拒绝，待下次确认）。
此前记录的「88 passed」是**半截同步**（缺 `filtering/geocoding/noise_filter/output/stc`
五个 `steps` 模块，但 `tests/test_filtering_kriging.py` 已上传 → 收集即失败）造成的
过期结论，已作废；远端因此长期不可用，直到 2026-10-29 19:20 整树重同步才恢复。
同场景 400×400×20 step3 实测：selection 0.429 s、network 0.0032 s、VCE 0.203 s、
estimate 0.0106 s、unwrap 1.562 s；ARM 单核性能明显低于本机 x86，但 640 核可并行拆分。
**2026-10-29 20:30 复测（新凭据、新代码树）**：pytest **122 passed + 1 skipped + 1 xfailed**
（skip = `reference/` 未同步的 `example.res` 用例；`GAP-S02` 当时尚未修复，故有 1 xfail），exit 0；
400×400×20 step3：selection 0.464 s、network 0.003 s、VCE 0.203 s、estimate 0.010 s、
unwrap 1.589 s（Npsc 520 / Narcs 1540 / 峰值 245.9 MB，与本地 w3 基准一致，无回归）；
`python -m depsipy.cli run --groups '[1:8]'` 在 60×60×15 合成工程上 rc=0，指标与本地相同
（sigma0 2.896、加密接受率 65.6%、`detrend_accepted=true`）。
同步命令：`~/.dsh_ssh/depsi-arm/cloud_test.sh`（已加 ControlMaster 单连接复用，同步 + pytest 一条命令）。

**接入方式（2026-10-29 20:00 前实测有效，后续凭据需重签）**：
跳板 `120.46.183.49:2234` 是 Go 写的 TCP 转发专用堡垒机——只允许 `direct-tcpip`、认证方法为
`none`，**凭据写在用户名里**：`-l 'jt_<id>:<token>'`。所以文档式
`ssh -J jt_id:token@120.46.183.49:2234 root@198.186.51.247` 不成立（OpenSSH 的 `-J` 不收密码，
会把整段当用户名却不带密码），必须用 `ProxyCommand` + `-W %h:%p` 并显式给 `-l`。
已把本机公钥 `~/.ssh/id_ed25519_depsi_arm` 装到目标 `root` 的 `authorized_keys`（目标侧免密，
公钥可随时撤销），辅助脚本在本机 `~/.dsh_ssh/depsi-arm/`（`run.sh` 执行、`push.sh` tar-over-ssh
同步——**本机无 rsync**、`jump.sh` 堡垒机转发、`jt_token` 单独存凭据、`pass.sh` 密码兜底）。
**已知不稳定 = 并发/会话数限制（用户 2026-10-29 19:41 判定；本轮到此停止，未再重试）**：
同一 token 开约 20 次会话后开始稳定返回 `platform authorization denied`（TCP 端口仍 OPEN、
SSH 握手与 `none` 认证都成功、只拒 `direct-tcpip`；冷却 ≥6 min 不复原）。
成因是配额、有效期还是并发限制**待确认**，但按并发限制处理。
**下次接入先做一件事**：给 `run.sh`/`push.sh` 加 `ControlMaster=auto` + `ControlPersist`
复用**单条**连接（一条连接跑多条命令），而不是每条命令重连——本轮就是被高频重连耗光的。
另外：目标密码 `123456`（用户 2026-10-29 临时设置）只在容器重建、公钥丢失时才需要。

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

9. 预检（`scripts/readiness.py preflight`）只读校验：参数必填/取值、清单 13 列、主影像唯一、
   裁剪交集、栅格字节数与声明尺寸、样点 NaN；返回 pass/blocked 与非零退出码。
   它发现了合成器的一个真缺陷：主影像曾同时写成干涉对与 `dummy` 行（`orbit+date` 重复），
   已修（现为 Nslc-1 辅影像 + 1 主影像，`truth.npz` 增加 `slave_index`）。
10. 原生 MATLAB 采集必须走 `matlab -batch` 且入口为函数 `depsi(param_file)`；`capture_baseline.m`
   检测到 Octave 直接报错，并在运行前剔除 `benchmarks/octave/shim`，避免 shim 顶替真实行为。

## 5. 分支规划

- 远程：`origin = https://github.com/Zoneip/depsi`（公开仓库，2026-10-29 创建）。
- `master`：当前 Python 实现（主开发线）。
- `rust`：预留分支，未来 Rust 版本在此线开发；创建于 2026-10-29，与 master 保持同一提交。
- 语言无关的规格（`specs/`）、规则（`rules/`）、对拍记录（`docs/parity.md`）与基准脚本
  保留在 master；Rust 分支按需同步，不单独复制一份。

## 6. 当前状态（每轮更新）

- 已完成·已验证：raw I/O、参数解析、27 模型、LAMBDA+ILS、估计几何、合成器、
  PSC/PSP 选择（精度 1.00 / 单元覆盖率 1.00，见 docs/selection-calibration.md）；
  当前测试集 124 项（本地 124 passed + 0 xfail；云端 122 passed + 2 skipped）
- 已完成·已验证（2026-10-29，M2 MVP）：`pipeline.py` + `cli.py`，`depsipy run <param_file>`
  串起处理组 2..8（selection → VCE/组网/估计/解缠 → 去趋势 → 加密 → 大气滤波 →
  噪声/STC → 结果落盘）。`tests/test_pipeline.py`（7 项）在 60×60×15 合成工程上验证：
  解缠相位 vs 真值中位误差 <1 rad、整数场 0 次修正、加密接受率 65.6% 且真 PS 富集、
  产品文件与 DePSI 列布局逐字节校验（`_psc_results_sel*` / `_ps_results_atmo_sel*`）。
  去趋势后若二次组网退化（整数修正达上限一半，或 sigma0 退化 >5×），主链自动回退到原始
  相位帧并记录 `detrend_accepted=false`（见 `docs/parity.md` P-11）。
- **Doris 产物直读 + 雷达地理编码（2026-10-29）**：用户确认 **Doris/RapidSAR 干涉代码可直接使用**
  （此前“许可不明、只参考”作废）。新增 `io/doris.py`（MFF 头、`.r00/.x00/.j00` 栅格、
  `.res` 元数据与精密轨道）、`core/orbit.py`（`intrp_orbit.m` 的多项式+样条状态）、
  `steps/geocoding.py` 的 `lph2xyz/xyz2t/xyz2lp/radar_to_geodetic`（`ps_output.m` 同款
  零多普勒+距离+椭球牛顿迭代）；`run_output` 在 `master_res` 非空时写 `_psc_latlon_sel*`
  并在 CSV 增加经纬度列。证据：`tests/test_doris_geocoding.py`（7 项，含真实
  `geocoding_v0.8/example.res` 解析与合成轨道往返 <1 mm）+ `test_pipeline` 集成。
- 已知未完成（**不得当作已支持**）：大气 Kriging 尚未把 PSC 的 APS 插值到 PSP（只有
  `steps/kriging.py` 的算法与单测）；加密只有单趟核心，缺幅度分箱迭代与 `mode_fix`；
  雷达地理编码**未在真实 `.res`+轨道上对拍**（只有合成轨道往返与真实 `.res` 解析）；
  状态落盘仍未实现（跨进程重启不可续算）。
- **已核实的缺口（strict xfail）**：无。三个缺口均已修复并有证据：
  `GAP-S01`（PSP 导出载荷，`test_psp_export_has_observation_payload`）、
  `GAP-U01`（孤立节点零列导致法方程奇异，`test_unwrapping_allows_isolated_node`）、
  `GAP-S02`（旁瓣去簇先取最弱 → 改为“最强为参考 + 连续相干邻域”；
  同数据 Octave 对拍 PSC/PSP 集合完全一致，见下条与 `docs/parity.md` P-15）。
- **MATLAB 运行环境盘点（2026-10-29）**：本机无 MATLAB（`command -v matlab` 为空、无安装目录），
  仅有 Octave 9.4；ARM 云端为 aarch64，官方 MATLAB Linux 版不支持（见上），故**当前唯一可跑
  MATLAB 原码的通道是 Octave + shim**。用户可提供 MATLAB 环境，接入前同数据对拍以 Octave 为准。
- **同数据 MATLAB(Octave) 对拍（2026-10-29）**：本机无 MATLAB、云端 ARM64 装不了官方 MATLAB；
  用 Octave 9.4 + shim 跑**未修改**的 `ps_selection.m`，与 depsipy 在同一 60×60×15 工程上比较：
  Npsc 139=139、Npsp 1456=1456、`(az,r)` 集合完全相同、相位最大差 2.5e-7 rad、h2ph 相同。
  工具：`scripts/octave_parity.py`；测试：`tests/test_matlab_parity.py`（无 Octave/工具箱时 skip）。
  对拍还查出并修复两个真 bug（`sidelobe_mask` 索引误用、`_decluster` 顺序反了）。
- 其他已核实缺口：`eachgrid`/`coherence_eachgrid` 选择分支、周期图估计、spider 的 Nparts 分区
  未实现或未对拍（详见 gap 矩阵 §2）。
- **上云前准备已完成（2026-10-29）**：`scripts/readiness.py`（environment/preflight/compare/
  from-mat/matlab）、`scripts/matlab/capture_baseline.m`、`scripts/matlab/export_comparison.m`；
  当时测试 107 passed + 3 strict xfail；现为 124 passed（无 xfail）。需求缺口见
  `docs/requirements-gap.md`，流程见 `docs/cloud-readiness.md`。
  **MATLAB 采集脚本未在真机执行**（无 x86-64 MATLAB 节点）。
- **性能对拍已完成（v0）**：estimate 987×（3000 弧×19 ifg，Octave 10.07s → 0.0102s，
  整数解 3000/3000 相同）、selection 10.8×（400×400×20，PSC 逐个相同）；详见
  `docs/performance.md` 与 `benchmarks/results/`。对拍还发现并修复 2 个真 bug（见 §4.8）。
- **向量化 W0-W2（step3）**：400×400×20 场景 selection 0.567→0.210 s、
  network 0.0065→0.0024 s、unwrap 1.332→0.957 s，step3 合计约 2.00→1.26 s；
  解缠 SuperLU `solve` 调用 29,298→1,542；详见 `benchmarks/results/w0-step3.json`、
  `w1-step3.json`、`w2-step3.json`。
- **W3（选择去簇的最近邻）已完成**：旁瓣去簇的最近邻改为增量维护 + k-d tree
  （大候选集）与阻塞精确回退（小集/平局兜底）。400×400×20 场景 selection 0.567→0.199 s
  （2.84×），**PSC 列表与旁瓣掩膜逐位相同**；step3 合计 2.008→1.285 s。
  详见 `benchmarks/results/w3-step3.json` 与 `docs/parity.md` P-10。
- **已回退的实验**：解缠校正循环的增量右端项（位相同但 1.09→1.66 s，构造稀疏列开销大于收益）。
  记录在 `docs/parity.md` P-10，代码已还原。
- 当前 step3 剩余瓶颈：`unwrap_network` 0.978 s（占 76%），其中
  `_resolve_integer_field_batched` 约 0.9 s，属 W 计划明确不向量化的顺序校正循环；
  若要继续压，需改算法（例如按列分块并行或减少校正轮数），需要另行评审。
- **W4 设计**：`specs/SPEC-003-稠密化向量化设计.md`（step6 稠密化，未实现生产代码）。
- **F1-F4 基础模块已落地（2026-10-29）**：`steps/densification.py`（KDTree 连接 + 批量 ILS +
  2of3/3of3 核心）、`filtering.py`、`kriging.py`、`noise_filter.py`、`stc.py`、`output.py`、
  `geocoding.py`；均有单测，**尚未接入 pipeline/端到端**。
- **三方 Python 复用分析**：`docs/python-reuse-matrix.md`（depsipy / 官方 Apache-2.0 /
  Doris 许可不明）；结论：depsipy 为主架构，官方只适配自包含 numpy/scipy 函数，
  Doris 仅参考、不复制代码；F1-F4 的模块级复用清单与 W 型迭代顺序已列。
- 授权：沿用「仅限科研，禁止商用」；2026-10-29 用户决策建立公开仓库（代码公开，许可不变）。
- 方案待评估：**MEX**（需 MATLAB + C++/OpenMP）与 **RunMat**（MATLAB 语法运行时，明确不支持 MEX）
  暂不下结论；MATLAB 环境接入后做真机基准与热点原型，再决定是否作为可选加速后端。

## 7. 纪律

- 改数值核前先跑 `pytest`；新增分支要有对应测试或不变量。
- 不与 MATLAB 逐位对齐为目标；对齐的是**数学模型**，差异写入 `docs/parity.md`。
- 不提交 `.venv/`、`*.raw`、`truth.npz`、`runs/`。
- 不擅自提交/推送/删除归档；需要时先问。

## 8. 待用户提供（阻塞项）

| # | 需要什么 | 用来决定什么 | 没它时不做的事 |
| --- | --- | --- | --- |
| 1 | 一份真实工程：`param_file.txt` + `input_file.txt` + 若干 `*.raw`（可只读） | 真实幅度/相位分布、`std_param` 是否合理、预处理链（RapidSAR）产物是否直接可读 | 调选择阈值、定性能常数、判断老版精度问题是否复现 |
| 2 | 服务器画像：未来提供 64 × 2 GB 显存大机器，算力不受限；**用户可 100% 提供 MATLAB 环境**。已实测一台 ARM64 Ascend 节点（640 vCPU / 1.9 TiB / 2×Ascend910），但官方 MATLAB 不支持 Linux ARM64；需 x86-64 节点才能装 MATLAB/MEX。MATLAB 版本与 Toolbox 待确认 | 并行策略、是否上 GPU/NPU 核、MEX 可行性 | 在拿到 x86-64 MATLAB 节点前不做 GPU/NPU 专用实现、不最终否决 MEX |
| 3 | 优先级：先"跑通出产物"还是先"补功能（加密/形变模型/大气分离）" | M2 与 M5 的顺序 | 继续 M2 之后的实现 |

已知可立即交付的东西（不需要新信息）：`depsipy` 的 raw I/O、参数解析、模型、LAMBDA/ILS、
选择、弧网、解缠、稠密化/滤波/噪声/STC/输出模块，以及可运行的 `depsipy run`（M2 MVP）。
