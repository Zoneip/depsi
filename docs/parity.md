# MATLAB ↔ depsipy 差异记录

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
