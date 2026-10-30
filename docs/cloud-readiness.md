# 上云前运行说明（MATLAB 基线 + Python 对照）

状态：**已准备，未在云端执行**（2026-10-29）。这里描述的命令都在本机通过了文件级与
失败路径测试；任何"已运行 MATLAB"的说法都必须以真机产物为准。

## 1. 本轮的目标与边界

只做三件事：① 记录环境；② 只读预检真实工程；③ 用**原生 MATLAB**采集基线并可与
Python 模块对照。

不做：补齐 Python 全流程、实现稠密化/大气/平滑算法、GPU/NPU 适配、大规模试跑。
这些要等生产入口、真实数据与 x86-64 MATLAB 节点到位后再定。

## 2. 已就绪的东西

```text
scripts/readiness.py               复查入口（environment / preflight / compare /
                                   from-mat / matlab）
scripts/matlab/capture_baseline.m  原生 MATLAB 采集（版本、Toolbox、函数解析、
                                   调用计数、产物清单、可选 profile）
scripts/matlab/export_comparison.m 显式导出对照数组（ids/epochs/fields）
tests/test_readiness.py            15 项：预检、对照、MATLAB 计划、失败路径
tests/test_interchange.py          2 项：MAT 转换、全 NaN 不算通过
tests/test_known_gaps.py           3 项 strict xfail：PSP 载荷/去簇/孤立节点
```

## 3. 标准流程

```bash
cd depsi
.venv/bin/python -m pytest tests -q            # 期望 85 passed, 3 xfailed

# ① 环境
.venv/bin/python scripts/readiness.py environment --output runs/env-$(date +%s).json

# ② 只读预检（不改动任何输入；--map 用于把生产绝对路径映射到本机副本）
.venv/bin/python scripts/readiness.py preflight /path/param_file.txt \
    --map pathmap.json --sample-pixels 256 --output runs/preflight-1.json

# ③ 生成 MATLAB 基线计划（默认不执行）
.venv/bin/python scripts/readiness.py matlab \
    --entry /path/to/ps_toolbox/main/depsi.m \
    --params /path/param_file.txt \
    --run-dir runs/baseline-1 --profile

# ④ 确认计划无误后再显式执行
.venv/bin/python scripts/readiness.py matlab ... --execute --timeout 7200
```

### 退出码

| 码 | 含义 |
| --- | --- |
| 0 | 通过 / 计划已生成 / 采集完成 |
| 1 | 对照不一致，或 MATLAB 执行失败 |
| 2 | 阻塞（缺文件、尺寸不符、参数非法、需要生产入口才能继续） |

`prepared` **不等于**跑过：`run.json` 里 `executed=false`、`status=prepared`。
只有 MATLAB 返回 0 且写出 `matlab-result.json` 且其中 `status=completed`，才会记 `completed`。
即便 `completed`，也只是"入口跑完"；是否真的执行了稠密化/滤波，要看 `profile` 采到的
函数调用计数（`report.calls`）。

## 4. MATLAB 侧的关键约束

1. **必须原生 MATLAB**。`capture_baseline.m` 检测到 Octave 直接报错；运行前会把
   `benchmarks/octave/shim` 从搜索路径里剔除，避免 `gcp/parpool/nanmean` 被替换。
2. **必须是 `matlab -batch`**，不接受交互式 session（无法保证复现）。
3. **入口必须是函数**，签名为 `depsi(param_file)`。若生产入口是硬编码 `work_dir` 的脚本，
   需要单独评审一个适配层，不能直接改生产文件。
4. **输出重定向**：计划会新建 `runs/<name>/work`，把输入清单第 6 列（`out` 前缀）改写到该目录，
   并把参数文件的 `input_file` 指向 `input.staged.txt`。原始参数与清单只读复制，绝不修改。
5. **相对路径**：清单里的相对栅格名会按"清单所在目录"解析成绝对路径并记录在
   `config.resolved_relative_paths`（DePSI 原本按工作目录解析，二者不同处会被显式列出）。
6. **第 9 步/第 6 步编号未定**时，`preflight` 会给出 `step_mapping_pending`；
   基线启动器直接判定 blocked，不会"猜"一个运行方式。

## 5. Python 对照怎么做

Python 侧目前只有选择、弧网、估计、解缠。对照按**模块**而非整链进行：

```bash
# MATLAB 导出（在采集脚本/适配层里调用，显式字段）
% export_comparison('sel_psc.mat', ids, epochs, struct('phase', P, 'amp_disp', D))
# 转成 NPZ，并写 sidecar JSON（units / reference / phase_direction）
.venv/bin/python scripts/readiness.py from-mat sel_psc.mat sel_psc.npz \
    --field phase --field amp_disp

# 对照
.venv/bin/python scripts/readiness.py compare mat.npz py.npz \
    --meta-left mat.json --meta-right py.json \
    --field phase --field amp_disp --atol 1e-9 --rtol 1e-7
```

规则：

- 先按**身份**对齐（PS 用坐标、弧用端点、时序用日期），再比数值。**点数相等不算通过**。
- 两侧 `units` / `reference` / `phase_direction` 必须一致，否则直接 `blocked`。
- NaN/Inf 掩膜必须逐位相同；全部为 NaN 的两侧判 `blocked`（没有可比样本），不算通过。
- 容差写进报告；尚未冻结业务容差，先用现有健康样本的误差量级（`1e-9` 相对量级）。
- 不因为"文件存在"就算实现：见 `docs/requirements-gap.md` 的功能矩阵（`◐` 条目仍未对拍）。

## 6. 服务器现状与限制

| 环境 | 事实 | 影响 |
| --- | --- | --- |
| 本机 | x86-64、16 vCPU、14 GiB、无 MATLAB | 可跑 Python 与 Octave 对拍 |
| ARM 节点 | aarch64（Kunpeng，640 vCPU / 1.9 TiB）、2×Ascend910、CANN 9.0、Python 3.10 | **官方 MATLAB 不支持 Linux ARM64**，不能装 MATLAB/MEX |
| 需要的节点 | x86-64 + MATLAB（含 Parallel Computing Toolbox） | 只有它才能出"真 MATLAB 基线" |

在拿到 x86-64 MATLAB 节点之前，`matlab` 子命令只能产出计划；这是设计上的阻塞，
不是缺陷。

## 7. 已知未验证项（不要当成已完成）

1. `scripts/matlab/*.m` 从未在真实 MATLAB 中执行过（只有语法与流程设计）。
2. `capture_baseline.m` 的 `profile` 调用计数依赖 MATLAB `profile` 的字段名，未实测。
3. 官方 Python 版（`reference/python/`）未安装、未运行；只做了静态对读。
4. 真实工程的幅度/相位分布、`std_param` 合理性、RapidSAR 产物兼容性均未验证。
5. 稠密化、大气滤波、Kriging、噪声滤波、STC、输出在地理与数值上都未实现。
