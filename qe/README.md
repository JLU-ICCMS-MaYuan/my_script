# QE 计算工具 - 使用与开发指引

## 状态概览
- CLI 框架位于 `qe/cli.py`，支持 `--json` 配置合并（优先级：命令行 > JSON > 默认），当前各子命令仅打印配置，实际工作流尚未落地。
- 预期工作流：结构优化 → 自洽 → 声子/电子/超导分析。对应队列、输入模板需要按需补充后接入 CLI。

## 快速上手（当前用于参数校验）
- 安装：`pip install -e .`。
- 命令示例（仅输出配置，不会提交计算）：
```bash
qe relax -i structure.cif -w ./relax --ecutwfc 80 --kpoints 16 16 16 --pp-dir ~/pot/qe_pp
qe phonon -i relax.out -w ./phonon --qpoints 4 4 4 --split --pp-dir ~/pot/qe_pp
qe superconductivity -i relax.out -w ./sc --method McAD --mu-star 0.1
```
- JSON 示例：`config_example.json`，使用 `--json config_example.json`。

## 代码结构速览
- `cli.py`：argparse 子命令（relax/scf/phonon/electron/superconductivity），负责合并 JSON + CLI 参数并打印配置。
- `pipelines/base.py`：批量处理基类，基于 `scheduler` 的 `Task` 与 `MixedParallelExecutor` 构建/运行任务，生成失败清单；可复用到未来的 QE 批量计算。
- `scheduler/`：任务封装与并行执行器；`core/types.py` 提供 TaskStatus 枚举。
- `config/` 与 `config_example.json`：赝势、k点、能量截断等参数示例。

## 开发与扩展指引
- 串联真实计算：为每个子命令实现对应 Workflow（如 `RelaxWorkflow`），负责生成输入、提交/监控作业，再在 `cli.py` 的 `command_*` 中调用。
- 批量模式：可复用 `pipelines/base.py`，在 `define_steps` 中声明步骤顺序与参数，再由 `scheduler` 调度。
- 队列支持：在 Workflow 中实现 bash/slurm/pbs 提交与完成判据，必要时加上重试与日志解析。
- 输入模板：如需标准 QE 输入，可在 `pipelines/templates/` 增加模板并在 Workflow 中渲染。

## 使用提示
- 赝势目录通过 `--pp-dir` 传入；当前仓库不会自动下载。
- 当实现真实计算后，建议在工作目录写入中间文件与最终报告，便于断点续跑与问题定位。
