# MATLAB ↔ depsipy 差异记录

状态：**持续记录（最后同步 2026-10-29）**。已记录 P-01…P-09；P-05、P-08 仍有开放待确认项。
格式：每条一节，写清 **MATLAB 行为 / depsipy 行为 / 原因 / 影响 / 是否可还原**。
只记录已核实的事实；待确认的写「待确认」。

---

## P-01 · `ps_model_definitions` 的全局 `althyp_index` 自增

- **MATLAB**：`case 'model'` 里对 `althyp_index = althyp_index + 1`（模型 ≤25 时）。
  这是副作用：同一个模型被第二次调用时断点会移动一位。
- **depsipy**：`build_functional_model(..., breakpoint=...)` 把断点作为显式参数，
  1-based 语义与 MATLAB 一致（断点位于 `breakpoint-1` 与 `breakpoint` 之间），不做自增。
- **原因**：原文的自增会让「多模型备选」时断点漂移，属于实现缺陷而非设计意图
  （`depsi.m` 里 `breakpoint` 是用户输入的固定时刻）。
- **影响**：使用 `breakpoint`/`breakpoint2` 且 `ps_model` 给出多个模型时，结果可能不同；
  单一模型时完全一致。
- **还原**：不提供。若需研究原文行为，`reference/matlab/main/ps_model_definitions.m` 可读。

## P-02 · LAMBDA 搜索步数上限

- **MATLAB**：`max_backtrack = n^3`，超过即 `break_flag=1` → 回退到 bootstrap 候选（非最优）。
- **depsipy**：默认 `max(100000, 1000*n^3)`；`LambdaModel(Qahat, compat=True)` 可还原 `n^3`。
- **原因**：该上限是性能保护，不是算法定义；实测它会让部分弧给出非最优整数解
  （`test_ils_is_optimal_in_the_ambiguity_metric` 会失败）。
- **影响**：depsipy 的整数解更优；`varfac` 分布略有变化（更好）。

## P-03 · VCE 的 trace 计算

- **MATLAB**：`N(k,j)=trace(QP·Qy1_k·QP·Qy1_j)` 按定义显式算，`O(Nsig²·Nifgs³)`。
- **depsipy**：利用 `Qy1_k` 是秩 1（`2·e_k e_kᵀ`，model==1 的 k=0 为 `2·ones`）的闭式：
  `N[k,j]=4·QP[j,k]·QP[k,j]`，`l[k]=2·(QP·y)_k²`，`O(Nsig²)`。
- **原因**：代数等价化简（`trace(A e_k e_kᵀ A e_j e_jᵀ) = A[j,k]·A[k,j]`）。
- **影响**：结果相同（数值上更稳，避免大矩阵连乘）；见 `tests/test_steps_network.py`。

## P-04 · 解缠法方程的求解方式

- **MATLAB**：每个 ifg、每轮修正都重算 `inv(Qyy)`、`chol(N)`，逐弧循环。
- **depsipy**：`scipy.sparse` 建 `AᵀWA`，`factorized()` 一次、多右端项求解；
  整数修正改为逐 ifg 的向量化循环（同一修正规则：取残差最大弧，`-round(e)` 步进）。
- **影响**：数学等价（同为加权最小二乘 + 同规则修正）；MATLAB 原文里
  「同一弧反复 ±1 的 skip 逻辑」用 `skip` 计数复现。

## P-05 · 分离网络（separate networks）的处理

- **MATLAB**：发现网络奇异时，把连接关系写回 `_psc_sel*.raw` 并**重新组网、重新估计**（最多 5 轮）。
- **depsipy**：直接按连通分量分组，每个分量各取一个参考点（主分量用全局参考或最小 `varfac` 弧，
  其余分量取最靠近中心的点），组内独立解算。
- **原因**：原文的重新组网在数学上等价于按连通分量分别解算，但会把整条链再跑一遍。
- **影响**：分组结果应一致（同组内参考点规则一致）；省掉重复的组网/估计轮次。
- **待确认**：原文对「孤立点」写回 `psc_array(:,2)=0` 后是否影响后续加密阶段的索引，
  需在 M5 做加密时核对。

## P-06 · 选区的旁瓣掩膜

- **MATLAB**（本工作区 hhf 版）：`strong_cand` = 3×3 局部极大且幅度 ≥ 95 分位；
  在行/列方向按「与更强候选相干 > 0.95 则剔除弱者」去簇，再做一次最近邻去簇；
  最终掩膜 = `strong_cand` 去掉存活者（即被抑制的旁瓣），被抑制像素的 `amp_disp` 置 NaN。
  另有 `weak_cand` 分支参与 `r_r/c_r` 拼接但**未参与最终掩膜**（死代码）。原文还有
  `side_lobe_threshold = inf;` 覆盖了上面按方法设定的阈值（同样是死赋值）。
- **depsipy**：`steps/selection.sidelobe_mask()` 复现「局部极大 + 95 分位 + 行/列相干去簇 +
  最近邻去簇」，忽略上述两处死代码。
- **影响**：掩膜应与原文一致；若校准后不同，记录在 `selection-calibration.md`。

## P-07 · 选择阶段的 I/O 与统计

- **MATLAB**：按网格缓冲区把整块 slc/ifg/h2ph 读入，逐缓冲区循环；`mean/std` 用 `nanmean/nanstd`。
- **depsipy**：一次幅度流式扫描累计 `sum`/`sum²`（等价于 `std` 的无偏式 `n/(n-1)` 修正），
  相位与 `h2ph` 只读取候选像素。
- **影响**：数值上等价（浮点求和顺序不同，差异 ~1e-7 相对）；I/O 从「全栈 × 缓冲次数」降为
  「全栈一次 + 候选数 × 景数」。分块边界不再影响结果（原文的旁瓣去簇按缓冲区进行，
  结果**依赖** `max_mem_buffer`，depsipy 不再依赖）。

## P-08 · `psc_distribution='uniform'`

- **MATLAB**：在平移半个网格的两套网格上各保留一个最优 PSC，并把参考点强制保留。
- **depsipy**：同样两套平移网格 + 强制保留参考点，但实现为两轮「每格保留最小值」的显式循环。
- **待确认**：原文对 `psc_new` 的排序/截断细节（`max_Npsc` 重排）是否影响最终点序；
  depsipy 最后按 `(az, r)` 排序，与原文 `sortrows` 一致。

## P-09 · 病态法方程的兜底

- **MATLAB**：`inv(C'*invQy*C)` 直接求逆；`h2ph` 近似常数时 cond ~1e10，参数可达 1e15。
- **depsipy**：`core/estimation._inverse()` 先算条件数，> 1e12 或非正定时退化为 `pinv` 并
  `RuntimeWarning`；正常情况用 Cholesky 解。
- **影响**：病态情形下 depsipy 给出最小范数解而不是垃圾值；健康情形数值一致。

## P-10 · 向量化实现的内部平局/顺序规则

- **MATLAB**：`ps_selection.m` 的旁瓣去簇用 `sort`/`sortrows`；本工作区版本按幅度升序处理，
  平局顺序未定义。解缠逐 IFG 循环，每个 IFG 独立修正整数。
- **depsipy**：`_decluster` 的 sweep 用 `np.lexsort((cand_amp, groups))` 先按组、再按幅度
  升序，组内等价于原 `argsort`，但平局时按原索引稳定排序；迭代最近邻 while 的裁决顺序不变
  （GAP-S02 已于 2026-10-29 修复，见 P-15）。解缠把外层 IFG 循环改为多右端批量推进，每个 IFG 的 skip/prev 规则
  不变。
- **影响**：只影响完全平局（幅度相同/残差相同）时的选点顺序，不改变数学模型；整数结果
  逐位一致，浮点相位按 `allclose` 验证。平局差异不提供还原开关。


## P-10 · 向量化后的平局与顺序规则（W1–W3）

纯 Python 循环改写为数组运算后，**平局由"输入顺序优先"统一决定**，并与参考循环逐位对齐：

| 位置 | 规则 | 验证 |
| --- | --- | --- |
| `_nearest_positive`（旁瓣去簇的最近邻） | 距离相同取 `pool` 中位置最小者；`d2 > 0`（坐标完全重合的候选**不**互为邻居），与参考 `argmin` 一致 | `test_decluster_matches_loop`（3 组 seed）+ 重复坐标用例 |
| `_decluster`（2026-10-29 重写） | 按 `ps_selection.m` 的 GexYuu 规则：每轮取**最强**候选为参考，剔除其两侧连续的相干邻域；不再按幅度升序 | `test_declustering_keeps_stronger_candidate`、`test_vectorization_parity.py`（3 项）、同数据 Octave 对拍（P-15） |
| `spider_arcs` | 弧对 `(min,max)` 去重后按行排序 | `test_spider_arcs_matches_loop` |
| `_resolve_integer_field_batched` | 每个 IFG 保留各自的 `skip/prev` 状态；批次只合并外层 | `test_batched_integer_correction_matches_loop` |
| `_solve_multi` | 分块只影响右端项切分，不改变解 | `test_solve_multi_matches_single_column_and_chunks` |

已知差异：`np.argsort` 不是稳定排序，因此凡依赖排序的步骤都改用 `np.lexsort`（显式给出
次序键）或显式 `argmin`（取首个最小值）。

**已尝试但回退**：把解缠校正循环的右端项改为增量更新（`rhs += AtW[:, idx] * step`）。
结果位相同，但 400×400×20 场景从 1.09 s 退化到 1.66 s——每轮构造稀疏列与 `sp.diags` 的
开销大于省下的稀疏乘；已回退，代码中不再保留。


## P-11 · 去趋势的平面估计保护

- **MATLAB**：`ps_detrend.m` 每次剔除残差超 100° 的一半点，直到剩余点数不足 3；在粗糙地形上
  会收敛到 3 个几乎共线的点，平面斜率可达真实值的 100 倍，随后 group4 的组网在失去
  模型一致性时失败。
- **depsipy**：保留"每次剔除超阈残差的一半"的规则，但**只要剩余点数少于初始点的 50%
  （或少于 10 个）就停止剔除**，用最后一次良态拟合。
- **影响**：平坦/常规地形下与 MATLAB 一致；粗糙地形下避免退化平面。`tests/test_pipeline.py`
  的 `test_detrend_keeps_the_phase_consistent` 断言平面系数 < 20 rad/px 且去趋势后网络
  0 次整数修正。
- **额外保护**：平面是按 IFG 估计的，并不保证与 `h2ph` 的缩放一致；在地形主导的数据上，
  二次组网可能因此失去整数一致性。`run_pipeline` 在去趋势后检查网络：若整数修正次数
  达到上限的一半，或 `sigma0` 退化超过 5 倍，则**拒绝去趋势结果并回退到原始相位**
  （`summary["detrend_accepted"]=false`）。

## P-12 · 解缠：孤立节点与删弧同步

- **MATLAB**：`ps_spatial_unwrapping.m` 先删掉连接数 ≤2 的弧（迭代），再删孤立点；法方程只含
  参考点固定的自由节点。
- **depsipy**：`_normal_matrix` 现在接收**自由节点掩码**（参考点 + 无弧节点都去掉），
  修复了孤立节点的零列导致 `Factor is exactly singular`；`remove_untestable` 删弧后，
  `dphase/dacheck/varfac` 用**删弧前**的弧键生成掩膜统一过滤，修复数组不同步（曾触发
  `IndexError`）。
- **验证**：`tests/test_known_gaps.py::test_unwrapping_allows_isolated_node`（原 strict xfail，
  现已转正）与 `tests/test_vectorization_parity.py`。

## P-13 · STC 的取数与单位

- **MATLAB**：`ps_spatio_temporal_consistency.m` 从 `_ps_results_*` 的 **defo（米）** 块读
  时间序列，`std` 直接用米，不做单位换算。
- **depsipy**：`pipeline.run_stc` 传入 `net.psc.defo`（米）并 `to_mm=False`；邻域用 KD-tree
  按 `(min_dist, max_dist)` 环带筛选，与 MATLAB 的网格 buffer 等价；邻居中含 NaN 的序列
  不参与投票（MATLAB 的 `_ps_results_*` 只含有效点，Python 侧把被删点保留为 NaN 行）。
- **影响**：STC 数值单位与 MATLAB 一致（米）；邻域实现方式不同但集合等价。


## P-14 · Doris 产物读取与雷达地理编码

- **Doris/RapidSAR**：MFF 文本头（`IMAGE_LINES`/`LINE_SAMPLES`）+ `.r00/.x00/.j00` 像素交错
  raw（与 DePSI raw 布局逐字节相同）+ `.res` 结果文件（`key: value`，含精密轨道
  `NUMBER_OF_DATAPOINTS` 与 `[time x y z]` 行）；读取脚本依赖 GDAL 的 MFF 驱动。
- **depsipy**：`io/doris.py` 自行解析 MFF 与 `.res`，栅格数据复用 `io/raw.py`，**不引入
  GDAL**；`core/orbit.py` 按 `intrp_orbit.m` 做 2 阶多项式（归一化时间）+ 解析速度/加速度，
  再用 not-a-knot 三次样条（MATLAB `interp1(...,'spline')`）插值状态；`steps/geocoding.py`
  按 `lph2xyz.m` 的零多普勒 + 距离 + 椭球三方程牛顿迭代（生产参数 MAXITER=100、
  CRITERPOS=1e-8），反方向按 `xyz2t.m`（迭代内线性插值、末次位置用样条）。
- **已知差异**：`metadataNew.m` 的 `last_pixel` 乘的是方位向多视因子 `container(8)`
  （`metadata.m` 用距离向 `container(9)`）；depsipy 采用 `metadata.m` 的正确因子，并在
  `docs/requirements-gap.md` 记录；`xyz2ell.m` 的闭式公式由迭代 `xyz_to_latlonh` 替代，
  差异 <1e-9 m。
- **验证**：`tests/test_doris_geocoding.py`（真实 `example.res` 解析、MFF 往返、
  合成直线轨道 `xyz2lp → lph2xyz` 往返 <1e-3 m、高度 500 m 往返）。


## P-15 · 同数据 Octave（MATLAB 原码）对拍

- **背景**：本机无 MATLAB；云端为 ARM64，官方 MATLAB 仅支持 x86-64。可运行 MATLAB 语法
  原码的唯一通道是 **Octave 9.4 + `benchmarks/octave/shim`**（7 个 shim；`ps_set_globals`
  的同名参数遮蔽问题已在驱动里显式赋值）。
- **工具**：`scripts/octave_parity.py` 用同一份合成工程生成 `sel_input.mat`（参数逐字段取自
  `param_file.txt`），分别跑 `ps_selection.m`（Octave）与 `depsipy.steps.selection`，比较
  PSC/PSP 计数、`(az,r)` 集合、相位与 `h2ph`；`tests/test_matlab_parity.py` 固化该对拍
  （无 Octave/无 `workspace/` 工具箱时自动 skip）。
- **结果（60×60×15，seed=2）**：Npsc 139 = 139、Npsp 1456 = 1456；PSC 与 PSP 的 `(az,r)`
  集合**完全相同**；相位最大差 `2.5e-7 rad`（float32 栅格存储），`h2ph` 完全相同。
- **对拍暴露并修复的真 bug**：
  1. `sidelobe_mask` 把 `_decluster` 返回的“保留索引”误当布尔掩膜（`suppressed[~keep]`），
     整数按位取反导致误杀/漏杀；PSP 曾因此少 53 个（1456→1403）。
  2. `_decluster` 的扫描顺序与参考相反（GAP-S02）：先取最弱而非最强，会连带剔除强点；
     已按 `ps_selection.m` 的“最强为参考 + 连续相干邻域”重写，并向量化（selection 400×400×20
     仍为 0.204 s，与 W3 的 0.199 s 持平）。
