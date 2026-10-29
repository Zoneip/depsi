# depsipy

DePSI v1.7.8.8（TU Delft Persistent Scatterer InSAR 时序列处理）的 Python 重实现。

目标：沿用同一套 `param_file.txt` / `input_file.txt` / `*.raw` 输入与文件格式，
跑同一条处理链，但显著更快、内存可控、可验证，并补齐原 Python 侧缺失的时序反演能力。

## 当前状态

| 模块 | 状态 |
| --- | --- |
| raw 栅格 I/O、参数解析 | 已完成 |
| 27 个功能模型、LAMBDA/ILS、估计几何 | 已完成 |
| 合成数据生成器、PSC/PSP 选择 | 已完成 |
| 弧网 + 解缠 | 代码完成，部分验证 |
| pipeline / CLI / 状态落盘 | 进行中 |
| 大气分离、加密、形变模型、噪声滤波、输出可视化 | 未开始 |

测试：`.venv/bin/python -m pytest tests -q`。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[fast,plot]'
```

可选依赖：`fast`（numba 编译核）、`plot`（matplotlib）、`hdf`（h5py）。

## 使用

```python
from depsipy.io import load_params, open_raw
from depsipy.io.params import parse_param_file
from depsipy.simulate import generate_project

# 读取 DePSI 工程参数
params = parse_param_file("param_file.txt")

# 以 memmap 方式打开无头复数 raw 栅格
stack = open_raw("slc.raw")

# 生成一个带真值的合成工程（用于验证）
project = generate_project(out_dir="/tmp/dp_sim", seed=0)
```

## 目录

```text
src/depsipy/
  io/        raw 栅格（memmap、1-based 裁剪、复数交错）、参数文件、工程状态
  core/      功能模型、LAMBDA/ILS、估计几何、编译核
  steps/     处理组：selection / network / unwrapping / ...
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

## 许可

仅限科研使用，禁止商用。算法与文件格式沿用 DePSI（Delft University of Technology）。
详见 [LICENSE](LICENSE)。
