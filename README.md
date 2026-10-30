# depsipy

DePSI v1.7.8.8（TU Delft Persistent Scatterer InSAR 时序列处理）的 Python 重实现。

目标：沿用同一套 `param_file.txt` / `input_file.txt` / `*.raw` 输入与文件格式，
跑同一条处理链，但显著更快、内存可控、可验证，并补齐原 Python 侧缺失的时序反演能力。

## 当前状态（2026-10-29）

| 模块 | 状态 |
| --- | --- |
| raw 栅格 I/O、参数解析、Doris 产物读取（MFF/`.res`/精密轨道） | 已完成 |
| 27 个功能模型、LAMBDA/ILS、估计几何 | 已完成 |
| 合成数据生成器、PSC/PSP 选择 | 已完成（M1 验收精度 1.00 / 覆盖率 1.00） |
| 弧网 + 解缠 | 代码完成，核心已验证（孤立节点等分支未覆盖） |
| 稠密化 / 大气滤波+Kriging / 噪声滤波 / STC / 输出+雷达地理编码 | **核心模块已落地，单测覆盖；除 Kriging→PSP 外已接入 pipeline** |
| pipeline / CLI | **MVP 已完成**（`depsipy run`，处理组 2..8，产品按 DePSI 布局落盘） |
| 状态落盘（跨进程续算） | 未开始 |

测试：`.venv/bin/python -m pytest tests -q` → 124 项全部通过（0 failed、0 xfail；
无 Octave/工具箱的环境会自动 skip 2 项：`reference` 示例与 MATLAB 同数据对拍）。

需求与缺口：[docs/requirements-gap.md](docs/requirements-gap.md)。
上云前流程：[docs/cloud-readiness.md](docs/cloud-readiness.md)。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[fast,plot]'
```

可选依赖：`fast`（numba 编译核）、`plot`（matplotlib）、`hdf`（h5py）。

## 使用

```python
from depsipy.io.params import load_params
from depsipy.io.raw import RawFile
from depsipy.io.project import read_input_file
from depsipy.simulate import generate_project

params = load_params("param_file.txt")          # 参数文件
stack = read_input_file(params.str("input_file"))  # 输入清单（已排序、含主影像）
with RawFile("slc.raw", stack.crop_final[1], "cpxfloat32") as f:
    block = f.read(1, 100, 1, 100)              # 1-based 裁剪读取

project = generate_project("/tmp/dp_sim", seed=0)   # 带真值的合成工程
```

完整处理链：

```bash
depsipy run param_file.txt --groups '[1:8]' --workdir runs/base1
```

`--groups` 接受 `3,4` 或 `[1:8]`；产品写在 `--workdir`（默认参数文件目录），文件名与列布局
沿用 DePSI（`_psc_results_sel*`、`_ps_results_atmo_sel*`、`_psc_atmo_filt_sel*` …）。
参数文件里给出 `master_res`（Doris `.res`）时，输出会额外写 `_psc_latlon_sel*` 并在 CSV
增加 `lat_deg/lon_deg/height_m` 列（零多普勒雷达地理编码）。
上云前的只读检查仍是 `python scripts/readiness.py --help`。

## 目录

```text
src/depsipy/
  io/        raw 栅格（memmap、1-based 裁剪、复数交错）、参数文件、工程状态
  core/      功能模型、LAMBDA/ILS、估计几何、编译核
  steps/     处理组：selection / network / unwrapping / densification / filtering / ...
  pipeline.py  端到端处理链（组 2..8）
  cli.py       `depsipy run`
  simulate.py  合成项目生成器
tests/       单元测试（pytest）
specs/       方案、ADR、验收标准
rules/       架构/编码/质量门规则
docs/        架构说明、对拍记录
benchmarks/  耗时与规模基准
```

依赖方向：`steps → core → io`。

## 文档

- [架构说明](docs/architecture.md)
- [与 MATLAB 原实现的差异](docs/parity.md)
- [性能对拍](docs/performance.md)
- [选择精度标定](docs/selection-calibration.md)
- [需求 ↔ 实现 ↔ 验证对照](docs/requirements-gap.md)
- [上云前运行说明](docs/cloud-readiness.md)

## 许可

仅限科研使用，禁止商用。算法与文件格式沿用 DePSI（Delft University of Technology）。
详见 [LICENSE](LICENSE)。
