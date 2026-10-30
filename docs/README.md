# docs/ · 索引

面向人的文档。规则在 `../rules/`，方案与 ADR 在 `../specs/`，归档在 `../achieve/`。
最后同步：2026-10-29（M2 MVP：`depsipy run` 端到端可用；reference 分为 matlab/python 两源）。

| 文档 | 内容 | 状态 |
| --- | --- | --- |
| [architecture.md](architecture.md) | 模块图、数据流、文件产物、状态格式 | 现行（状态落盘未实现） |
| [parity.md](parity.md) | MATLAB ↔ depsipy 的逐条差异与理由 | 已记录 P-01…P-13；P-05/P-08 有开放待确认项 |
| [selection-calibration.md](selection-calibration.md) | 合成数据上选择参数校准记录 | M1 验收通过（2026-10-29） |
| [performance.md](performance.md) | 四版本可测性、两场景实测、机制归因、口径与免责 | v0 实测（2026-10-29） |
| [requirements-gap.md](requirements-gap.md) | 需求原文、业务步骤↔源码调用、功能矩阵、参数待核 | **现行**（2026-10-29） |
| [python-reuse-matrix.md](python-reuse-matrix.md) | depsipy / 官方 Python / Doris 三方代码差异与复用清单 | **现行分析**（F1-F4 已落地并接主链） |
| [cloud-readiness.md](cloud-readiness.md) | 上云前流程：环境/预检/原生 MATLAB 基线/对照规则 | **现行**（未在云端执行） |

待建：`user-guide.md`（真实工程完整步骤）、`benchmarks/report.md`（M4 规模基准报告）。
