# 性能实测报告（v0，2026-10-29）

状态：**v0 实测完成**。本报告全部数据来自本机（16 vCPU / 14 GiB / 无 GPU）。
用户已告知未来提供 64 × 2 GB 显存的大机器（型号与加速栈待确认），届时需按本文 §7 复测；
用户另告知**可 100% 提供 MATLAB 环境**（版本/接入时间/Toolbox 待确认），届时可用真机 MATLAB
替换 Octave 代跑，并补 MEX 基准。在拿到机器与 MATLAB 前不做 GPU 专用实现、不最终否决 MEX。

> 口径：**实测**与**估算**严格分开。所有"实测"都可按本文末命令复现；凡涉及 MATLAB
> 的地方都写明是用 **Octave 9.4** 代跑（本机无 MATLAB），并给出外推区间。

## 0. 环境

| 项 | 值 |
| --- | --- |
| CPU / 内存 | 16 vCPU / 14 GiB，无 GPU |
| 系统 | Linux（Debian 13 内核源），`df` 可用 203 GB |
| MATLAB 版本 | **M 文件用 Octave 9.4.0 代跑**（含 image/statistics/parallel 包） |
| depsipy | Python 3.13.5、numpy 2.5.4、scipy 1.18.1、numba 0.68.1 |
| 数据 | `simulate.generate_project` 合成：DePSI 原始格式 + 真值 |
| 算力规划 | 未来 64 × 2 GB 显存机器（待接入，未参与本报告） |

## 1. 四个版本能不能测

| 版本 | 类型 | 能否测速 | 原因 |
| --- | --- | --- | --- |
| `ps_toolbox/ps_toolbox`（v1.7.7.5） | MATLAB | ✔ 用 Octave | m 代码可直接跑，缺的函数用 shim 补齐（见 §5） |
| `ps_toolbox_v1.7.8.8_hhf` | MATLAB | ✔ 用 Octave | 同上；其 `parfor` 在 Octave 中**串行执行** |
| `Doris干涉加时序depsi/` | Python 干涉链 | ✘ | 依赖 DORIS 可执行文件与 SAR 原始数据，本机两者皆无；且它不含时序反演 |
| `InSAR_affair-shallow-2026-10-19/` | Python DANI-NET | ✘ | 原工作区副本已于 2026-10-29 移除；本机也未装 torch/NPU 训练环境，且与 PS 时序列不是同一任务 |
| `depsipy`（本次新写） | Python | ✔ | — |

**结论**：能做同类比较的只有两个 MATLAB 版本与 depsipy。另外两个不是"同一程序的四个版本"，
而是两条独立产品线（干涉链 / 深度学习），把它们放进速度对比没有意义——这一点先说明，
避免得出错误结论。

## 2. 实测一：逐弧模糊度估计（LAMBDA / bootstrap）

场景：**3000 弧 × 19 干涉图**，同一份 `est_input.mat` 喂给三方（同一个 `dphase`、`h2ph`、
`Btemp`、`Bdop`、`std_param`，`ps_method='boot'`，`psc_model=[1]`）。

| 实现 | 耗时（热） | 吞吐 | 峰值内存 | 整数模糊度一致数 |
| --- | --- | --- | --- | --- |
| v1.7.7.5 @ Octave | 10.23 s | 293 弧/s | — | 基准 |
| v1.7.8.8 @ Octave | **9.82 s** | 305 弧/s | 112.6 MB | 与 v1.7.7.5 完全相同 |
| **depsipy（numba）** | **0.0100 s** | **300,000 弧/s** | 158 MB | **3000/3000** |
| depsipy（纯 Python 兜底） | 3.17 s（按 300 弧外推） | 946 弧/s | — | 300/300 |

- 相对两版 MATLAB：**约 980×**（numba）/ **3.2×**（不装 numba）。
- 数值一致性（3000 弧）：`acheck` 最大差 **0**；`varfac` max|Δ| = **1.8e-14**；
  `ens_coh` max|Δ| = **6.5e-15**；中位数完全相同（varfac 1.8196、coh 0.8199）。
- 冷启动：首次调用含 numba JIT，约 **0.19 s**，之后进程内复用（`cache=True` 也会落盘）。

> 两个 MATLAB 版本在估计上耗时相同、结果相同：它们只差"日志 + CPU 调度块"，
> 而调度块的 `parfor` 在 Octave 中串行。真机 MATLAB 用 14 worker 可能更快，
> **但这一点本机无法验证**。

## 3. 实测二：PSC / PSP 选择（含 I/O）

场景：**400×400 像素、20 景**（19 干涉图 + 主影像）、真 PS 占比 10%、
网格 18 m、`psc_threshold=0.3`、`psc_distribution='nonuniform'`、启用旁瓣掩膜。
真值来自生成器，可算"选出来的点里有多少是真 PS"。

| 实现 | 耗时 | PSC | PSP | 峰值内存 | PSC 精度 |
| --- | --- | --- | --- | --- | --- |
| v1.7.7.5 @ Octave | 13.45 s | 451 | 5,589 | — | **0.572** |
| v1.7.8.8 @ Octave | 5.52 s | 520 | 55,207 | 387.7 MB | **1.000** |
| **depsipy** | **0.53 s** | **520** | 52,671 | **82.5 MB** | **1.000** |

- 相对 v1.7.8.8：**10.4×** 更快、**4.7×** 更省内存，且**选出的 520 个 PSC 与它逐个相同**。
- 相对 v1.7.7.5：**25.4×**。
- 顺带结论：**老版不只慢，还会选错**——精度 0.572 vs 新版 1.000，差别来自 hhf 分支新增的
  旁瓣去簇（行/列相干 + 最近邻）。这是个实质性修复，不是风格差异。

## 4. 为什么快（机制，逐条可查）

| # | 机制 | 证据 |
| --- | --- | --- |
| 1 | 选择阶段只读候选像素，而非整块缓冲 | 同一 3000×3000 文件：整幅 memmap 读 32.8 ms，稀疏读 1% 像素 **4.6 ms（11×）**；`freadbk` 式逐行 seek+read 51.2 ms |
| 2 | 估计的几何与 LAMBDA 分解按**模型**算一次，逐弧复用 | `Geometry` / `LambdaModel`；MATLAB 每弧重算 `inv(Qy)`、`inv(C'invQy C)`、`decorrel` |
| 3 | LAMBDA 分支定界编译化 + 并行 | numba `prange`；单弧 3.3 µs |
| 4 | VCE 的 trace 用秩 1 闭式 | `trace(A e_k e_kᵀ A e_j e_jᵀ) = A[j,k]·A[k,j]`，O(Nsig²·Nifgs³) → O(Nsig²) |
| 5 | 解缠法方程因子化一次、多右端项 | `scipy.sparse.factorized`；MATLAB 逐弧/逐轮重算 `chol(N)` |
| 6 | 状态用 npz+json 而非 `save/load *_project.mat` + `pack` | 每组之间不再整工作区落盘 |

## 5. 口径与免责（必须一起读）

1. **Octave ≠ MATLAB**。m 代码循环上 Octave 通常比 MATLAB JIT 慢 **2–10×**。
   因此"980×"是**相对 Octave 的实测**；相对真机 MATLAB 的合理区间是 **100–500×**（估算，未实测）。
   判断依据：MATLAB 没有"把标量循环编译成机器码"的能力，量级差异来自解释执行 vs 编译。
2. **parfor**：MATLAB 侧可开多 worker；本机 Octave 串行，所以 MATLAB 的实测可能优于本表。
   若要给用户一个"最保守"的结论：**depsipy 估计步骤至少与开了 parfor 的 MATLAB 同量级或更快**，
   因为 depsipy 的 0.0100 s 已经包含了 16 线程。
3. 选择的 I/O 优势与解释器无关（是字节数），这部分结论对 MATLAB 同样成立。
4. 合成数据 ≠ 真实数据：真实 SLC 的幅度分布、失相干、地形相位都会改变耗时的常数项
   （但改变的是常数，不是复杂度）。
5. 本报告是 **v0**：只覆盖估计与选择两步；pipeline/CLI 尚未落地（M2），端到端数据待 M2 后再补。

## 6. 对拍发现的两个真 bug（已修，有回归测试）

用 Octave 跑 MATLAB 原码、与 depsipy 逐量对拍（`Qy1`/`Qahat`/`D`/`Chi2`/`varfac`/`ens_coh`）
才发现的两个问题：

| # | 症状 | 根因 | 影响 | 回归测试 |
| --- | --- | --- | --- | --- |
| 1 | `varfac` 中位 0.165 vs MATLAB 1.631（差 9.87×） | 度转弧度写错：`(pi*15/180)^2` 被实现成 `(pi*15*(pi/180))^2` | 方差分量整体放大 9.87×，权重与阈值判断随之偏移 | `test_default_variance_components_match_matlab` |
| 2 | 200 弧里 40 弧整数模糊度错 1 个；批处理结果呈线性斜坡 | numba `prange` 循环**外**声明临时数组 → 多线程共享互相覆盖 | 估计结果错乱且不稳定（同样输入不同结果） | `test_batch_resolution_matches_per_arc_reference` |

这两类问题"读代码看不出来、单测也未必抓到"，只有**对拍**能发现。因此建议把 Octave
保留为常设对拍通道（本机已装 9.4.0）。

## 7. 复现命令

```bash
cd depsi
export PYTHONPATH=src

# 1) 生成场景（同一份数据喂三方）
.venv/bin/python benchmarks/make_scenarios.py --out /tmp/dpbench \
        --lines 400 --pixels 400 --nslc 20 --arcs 3000

# 2) MATLAB 侧（Octave 代跑）
cd benchmarks/octave
octave --no-gui --quiet bench_estimation.m <toolbox_dir> /tmp/dpbench
octave --no-gui --quiet bench_selection.m  <toolbox_dir> /tmp/dpbench /tmp/selwork

# 3) depsipy 侧
cd ../..
.venv/bin/python benchmarks/run_bench.py --scenario /tmp/dpbench

# 4) 单弧对拍（打印 MATLAB 中间的 Qahat/D/Chi2/varfac/ens_coh）
cd benchmarks/octave
octave --no-gui --quiet dbg_core.m <toolbox_dir> /tmp/dpbench
```

`<toolbox_dir>` 取 `workspace/ps_toolbox_v1.7.8.8_hhf` 或 `workspace/ps_toolbox/ps_toolbox`。

## 8. Octave 侧需要的 shim（记录）

MATLAB 有、Octave 没有的函数，放在 `benchmarks/octave/shim/`：
`gcp`、`parpool`、`feature`（`feature('numCores')`）、`nanmean`、`nanstd`、`nanmax`、`nanmin`。
另外 `ps_set_globals` 会把同名**参数**遮蔽成空全局变量，所以驱动脚本必须把
`Npsc_selections`、`weighting`、`varfac_threshold` 等全部显式赋值（`depsi.m` 里由
`ps_readinput_parameters` 完成）。这两点是跑通 MATLAB 码的关键，已写进脚本注释。
