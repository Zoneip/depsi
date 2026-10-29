# 选择步骤校准记录（M1）

数据：`simulate.generate_project(Nlines=120, Npixels=120, Nslc=12, seed=0, ps_fraction≈10%,
noise_amp=0.25, DS 幅度调制 σ=0.45, PS 幅度调制 σ=0.12)`，参数取自生成器写出的 `param_file.txt`
（`psc_threshold=0.3`、`psc_selection_gridsize=18 m`、`psc_distribution='uniform'`、
`do_apriori_sidelobe_mask='yes'`、`peak_tolerance=0.7`）。

## 结果

| 指标 | 实测 | 说明 |
| --- | --- | --- |
| PSC 数 | 48 | 网格 12×4，一格一个，与 `psc_distribution` 语义一致 |
| **精度**（选出的是真 PS） | **1.00（48/48）** | 假阳性 0 |
| **单元覆盖率** | **1.00** | 含「可达真 PS」的网格单元全部产出 PSC |
| 幅度离散（真 PS） | 中位 **0.117**，90 分位 0.13 | 与合成注入的 σ=0.12 一致 |
| 幅度离散（真 DS） | 中位 **0.426** | 与注入的 σ=0.45 一致 |
| 旁瓣掩膜抑制像素 | 702，其中真 PS 占比 **0.00** | 未误杀 PS |
| 选区耗时 | 0.03 s（120×120、12 景） | 流式幅度 + 稀疏相位读取 |

## 结论与对原实现的解释

1. **召回率（像素级）不是合适的指标**：DePSI 按网格把点数压到「每格一个」，
   像素级召回天然是 `格数 / 可达 PS 数`。本轮 1514 个可达真 PS 只保留 48 个，
   这正是网格化的设计意图。因此验收改用 **单元覆盖率**（≥95%）+ **精度**（≥95%）。
2. `psc_threshold=0.3` 在合成数据上分离良好（PS 中位 0.117 vs DS 0.426），
   阈值有 ~2.5× 余量；真实数据若 PS 幅度波动更大，应先用历史工程调阈值。
3. 旁瓣掩膜（局部极大 + 95 分位 + 行/列相干去簇 + 最近邻去簇）抑制了 702 个像素，
   均为强 DS 的旁瓣，未触及真 PS。参数 `peak_tolerance` 从原代码的 0.70 沿用。
4. 原文两处死代码（`weak_cand` 未参与掩膜、`side_lobe_threshold = inf`）已按
   `docs/parity.md` P-06 忽略。

## 复现

```bash
cd depsi
PYTHONPATH=src .venv/bin/python -m pytest tests/test_simulation_selection.py -q
```

## 待办

- 真实工程（用户 Spacety 数据）上的阈值复核，需要 `psc_threshold` 与真实幅度离散分布；
  缺少数据前不做进一步调参。
