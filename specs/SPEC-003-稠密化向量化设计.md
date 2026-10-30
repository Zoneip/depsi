# SPEC-003 · 稠密化（step6）向量化设计

状态：**设计稿（2026-10-29，W4 交付；未实现生产代码）**
关联：`workspace/Depsi效率提升需求.docx` R1（step3/step6 效率低）、`SPEC-001`、本计划
`../.trae/documents/for循环向量化-W型迭代计划.md`

## 1. 背景与现状

需求清单明确指出 step3 与 step6（稠密化）是主要低效步骤。MATLAB
`ps_densification.m`（1282 行）有 21 个 `for` + 1 个 `while`，densification 在
group5（atmo 分支）和 group6（defo 分支）各调用一次；当前 `depsipy` 尚未实现该模块。

本设计只定义**向量化优先的实现骨架**，不改变 `2of3/3of3/mode_fix` 等判定语义；
生产代码在 M5 落地。

## 2. 数据模型（实现时必须先建立）

以一次 buffer / 一个 selection 为单位，所有候选与 PSC 用数组表示，禁止
在 `for` 内反复用 `find`/`ismember` 取子集：

| 名称 | 形状 | 说明 |
| --- | --- | --- |
| `psp_az`, `psp_r` | `(Npsp,)` | PSP 像素坐标（1-based） |
| `psp_phase`, `psp_h2ph` | `(Npsp, Nifgs)` | 候选相位与 h2ph 因子 |
| `psp_amp_disp` | `(Npsp,)` | 幅度离散 |
| `grid_id` | `(Npsp,)` | 预计算 grid cell id（az cell × Ngrid_r + r cell） |
| `ref_candidates` | `(Ngrid, Nest)` | 每个 grid 最近的 `Nest` 个 PSC 索引（cKDTree 一次查询） |
| `ps_az`, `ps_r`, `ps_phase_unw`, `ps_acheck` | `(Nps,)` / `(Nps,Nifgs)` | 当前已加密 PS |
| `amp_disp_bins` | `(Ndens_iterations,)` | 每轮幅度分箱阈值 |

## 3. 向量化流水线

### 3.1 候选提取与 grid 映射

MATLAB `L334/L400/L982` 的 buffer/grid/迭代三重循环改为：

```python
mask = amp_disp_grid < amp_disp_bins[p]
flat_idx = np.nonzero(mask)[0]
cell_a = np.searchsorted(grid_az_starts, rows[flat_idx] + 1, side="right") - 1
cell_r = np.searchsorted(grid_r_starts, cols[flat_idx] + 1, side="right") - 1
cell_id = cell_a * Ngrid_r + cell_r
```

一次得到所有候选的 cell id；后续按 `cell_id` 分组，不再逐 cell `find`。

### 3.2 每个 grid 的参考 PSC

MATLAB `L416-L420` 对每个 grid 中心排序全部 PSC。改为：

```python
centers = np.column_stack([center_az * az_spacing, center_r * r_spacing])
ref_dists, ref_idx = cKDTree(psc_xy).query(centers, k=Nest)
```

`ref_idx` 形状 `(Ngrid, Nest)`；若 `Nest > Npsc`，用固定哨兵并在估计阶段掩膜。

### 3.3 候选参数估计

MATLAB `L459-L483` 对每个 grid、每个候选、每个参考点调用一次 ILS/perio。改为：

```python
# 对每个候选 × Nest 展开为弧： (Npsp * Nest, Nifgs)
dh2ph = 0.5 * (psp_h2ph[cand][:, None, :] + psc_h2ph[ref_idx[cell_id]])
dphase = psp_phase[cand][:, None, :] - psc_phase[ref_idx[cell_id]]
res = estimate_arcs(dh2ph.reshape(-1, Nifgs), ..., models=ps_model)
```

`estimate_arcs` 已是 numba 批量核，一次处理全部候选，消除逐 grid/逐候选的
Python 调用与 MATLAB 标量循环。

### 3.4 `2of3` / `3of3` / `mode_fix` 判定

MATLAB `L488-L552` 逐 PS 设 `error_flag`。改为数组位掩码：

```python
diff12 = np.abs(xcheck[:, :, 1] - xcheck[:, :, 0]).sum(axis=1)
diff23 = np.abs(xcheck[:, :, 2] - xcheck[:, :, 1]).sum(axis=1)
diff31 = np.abs(xcheck[:, :, 0] - xcheck[:, :, 2]).sum(axis=1)
flag = (diff12 != 0) * 1 + (diff23 != 0) * 2 + (diff31 != 0) * 4
# flag == 0/3/5/6/7 等分支用 np.where / 掩膜索引
```

修正 `psp_acheck/psp_param` 用 `flag` 掩膜批量赋值，不逐 PS 循环。
`mode_fix` 的众数/修复规则单独写成一张掩膜表，先与 MATLAB 小样本对拍。

### 3.5 最终 PS 拟合与 VCE

MATLAB `L638-L697` 逐 PS 调 `rhs*y`。改为：

```python
xhat = rhs @ Y.T              # (Npar, Nps)
yhat = B1 @ xhat
ehat = Y.T - yhat
ens_coh = np.abs(np.exp(1j * ehat).mean(axis=0))
```

VCE 分支的 `sig2hat`、`Qahat` 也用矩阵运算；`Qahat` 只对需要输出的 PS 计算，
或分块计算，避免 `Nps×Npar²` 全量存储。

### 3.6 最近邻加密

MATLAB `L1025-L1028` 对每个候选求到已有点的最小距离。改为：

```python
_, nearest = cKDTree(ps_xy).query(psp_xy, k=1)
```

一次查询；`psp_dh2ph_grid`、`psp_dphase_grid`、`ps_acheck_ref` 用
`nearest` 数组索引构造。

### 3.7 外层迭代

`Ndens_iterations` 的 bin 循环保留（语义是按阈值逐轮放开候选）；
`sweep` 内所有逐 grid/逐候选/逐 PS 循环数组化。`while model_index<=Nmodels`
保留（模型选择顺序依赖）。

## 4. 预期函数与文件

M5 实现时建议新增 `src/depsipy/steps/densification.py`：

- `grid_ids(rows, cols, grid_az, grid_r)`
- `nearest_ps(ps_xy, psp_xy, k)`
- `build_densification_arcs(psp, ref_idx, cell_id, ...)`
- `densify(psd, options) -> DensificationResult`
- `check_consistency(xcheck, mode)`（`2of3`/`3of3`/`mode_fix`）

并复用：
- `core.estimation.estimate_arcs`（批量 ILS）
- `steps.network.form_network`（加密后网络）
- `steps.unwrapping.unwrap_network`（加密后解缠）

## 5. 验收（实现时）

1. 小数据（100×100、12 景、`Ndens_iterations=2`、`Nest=3`）：
   - PS 数、`acheck`、`ps_param`、`ens_coh` 与 MATLAB/循环参考版一致；
   - `2of3`、`3of3`、`mode_fix` 三种分支各有一组用例。
2. 扩展数据（400×400、20 景）：
   - 候选数翻倍时，内层耗时接近 O(Npsp log Npsp)，而非 O(Npsp²)；
   - 峰值内存 ≤2× 循环参考版。
3. 回归：`pytest tests -q` 全绿；不改变现有 step3 结果。

## 6. 风险与未决

- `psp` 数量可能非常大，`Npsp×Nest×Nifgs` 展开数组需分块；按 buffer/grid 块限制内存。
- `ref_idx` 的 `Nest` 不足、重复坐标、NaN 相位需显式掩膜；MATLAB 的 `sort` 平局规则需对拍。
- `mode_fix` 分支语义复杂，必须先用 MATLAB 小样本锁定，不允许“顺手改”。
- group5 与 group6 的 densification 参数不同（atmo vs defo），设计需两份验收数据。
