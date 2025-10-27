# QE 模块化与解耦规划

## 总体目标
- 建立清晰分层的 QE 计算脚本体系：CLI/配置 → 工作流 orchestration → 输入/输出模板 → 作业提交 → 数据后处理。
- 将目前通过 `kwargs` 动态注入的属性改为强类型参数对象，提升健壮性与可读性。
- 按“结构优化 / 自洽 / 声子 / 电声耦合 / 电子性质 / 超导 / 批量工具”等场景划分模块，避免巨型类承担所有职责。
- 统一日志、错误处理、任务监控方式，并为后续绘图、数据导出提供稳定接口。

## 当前痛点
- `qe_inputpara.py`、`qe_writeinput.py`、`qe_writesubmit.py`、`qe_run.py` 体量巨大，职责交叉，`if-elif` 零散处理所有模式。
- 参数管理依赖 `hasattr` + `kwargs`，缺少默认值校验与类型提示；不同工作流的特有参数混杂在同一类中。
- 输入文件、提交脚本通过字符串拼接生成，模板重复多、难以维护；缺乏统一的输出目录/命名策略。
- 作业提交逻辑 (`qe_submitjob.py`) 同时覆盖 bash/slurm/pbs/lsf，分支嵌套，难以扩展到新的调度器或监控策略。
- 数据处理（如 lambda、phonon DOS、HS 路径）散落在 `qe_base` 与子类中，缺乏统一接口，难以复用到可视化或导出。
- CLI (`set_args.py`) 将每个子命令硬编码指向类，耦合度高，文档和帮助信息无法共享默认值/参数说明。

## 分层与目录建议
```
qe/
  core/
    config.py          # ArgumentParser 解耦，加载默认配置
    logging.py         # 统一日志初始化
    structures.py      # 结构读取、倒格矢、赝势选择
    parameters/
      base.py          # 通用 dataclass + 验证逻辑
      relax.py / scf.py / phonon.py / electron.py / superconduct.py / batch.py / sctk.py
  workflows/
    base.py            # 抽象 Workflow，提供 run() / prepare() / submit() 钩子
    relax.py, scf.py, phonon.py, ...
  io/
    writers/
      pw.py / ph.py / matdyn.py / job_script.py ...
    templates/         # 可选：Jinja2 or str.format 模板
    readers/
      phonon.py / electron.py / lambda.py
  submission/
    backends/
      bash.py / slurm.py / pbs.py / lsf.py
    monitor.py         # PID/JobID 监控、错误识别
  postprocess/
    phonon.py / electron.py / superconduct.py / visualization.py
  cli/
    __main__.py        # `qe_main.py` 精简为入口，依赖 core.config
    commands.py        # 子命令注册表
  utils/
    files.py / shell.py / validators.py
  legacy/
    eliashberg/        # Fortran 代码保持隔离，提供包装器
```

## 参数与配置管理
- 为通用参数定义 `BaseParameters` dataclass，负责工作路径解析、结构备份、赝势准备。
- 每个工作流继承 `BaseParameters` 并新增字段（如 `PhononParameters.tr2_ph`、`SuperconductParameters.screen_constants`），集中默认值与验证逻辑。
- 使用工厂函数/注册表将 CLI 子命令映射到 `Parameters` + `Workflow`，方便新增模式。
- 引入 `from_dict` / `from_args` 静态方法，替换当前 `init_from_config` + `kwargs` 方案。

## 工作流拆分
- 创建 `Workflow` 抽象类：定义 `prepare_inputs()`、`prepare_submission()`、`execute()`、`postprocess()` 等接口。
- `RelaxWorkflow` / `SCFWorkflow` / `PhononWorkflow` / `ElectronWorkflow` / `SuperconductWorkflow` / `BatchWorkflow` / `SCTKWorkflow` 分别实现上述接口，根据模式选择具体的 writer、submitter、postprocess 组件。
- 将当前 `qe_run.py` 中的模式分支改为注册机制：`WORKFLOW_REGISTRY["relax"] = RelaxWorkflow`，命令行只负责查表执行。

## 输入/输出生成
- 将 `qe_writeinput.py` 拆成按程序归类的 writer（如 `PwInputWriter`、`PhInputWriter`、`MatdynInputWriter`、`DosInputWriter`）。
- 采用模板驱动（Jinja2 或 Python f-string 模板），模板中只填充参数对象，避免重复拼接。
- 把 `ATOMIC_SPECIES`、`CELL_PARAMETERS` 等常见片段下沉为复用的 helper。
- 针对 `split_dyn0` 等批量场景，引入 `BatchWriter` 管理子目录创建和文件命名。

## 作业提交与监控
- 设计 `SubmissionBackend` 接口：`submit(job_script)` 返回 jobid；不同系统各自实现（bash/slurm/pbs/lsf）。`qe_submitjob.py` 简化为选择后端、统一错误处理。
- 将 `submit_mode0/1/2/3` 拆成“单 job / 持续监控 / 分片任务”策略类，组合 `SubmissionBackend` 使用。
- 统一监控工具：轮询 `jobid`、捕获输出日志，错误检查放入独立 `LogInspector`。`checksuffix`、`checkerror` 等方法移动至此。

## 数据处理与后处理
- 把 `qe_base` 中的结构处理、lambda 计算、HS 路径等方法拆成独立模块：
  - `core.structures.StructureAnalyzer`：结构、赝势、对称性、k/q 网格。
  - `postprocess.phonon.LambdaCalculator`：读 `gam.lines`、`freq`、`a2F`。
  - `postprocess.paths.HighSymmetryPath`：支持从 `matdyn.in`、ASE 自动推断或手动输入。
- 提供统一的数据导出接口：`to_dataframe()`、`to_csv()`、`to_plot_model()`，后续绘图/可视化共享。

## CLI 与用户交互
- 精简 `qe_main.py`：入口只负责解析参数、加载默认配置、初始化日志、调用对应 Workflow。
- `set_args.py` 拆成 `cli/commands.py`，通过装饰器或注册函数集中定义命令、描述、默认值；帮助信息自动从参数对象注释生成。
- 将 `input()` 交互（如选高对称路径）改为 CLI 参数或独立交互命令，避免在核心逻辑里阻塞。

## 工具与公共功能
- 将 `qe_base` 拆为多个 `utils` 模块（文件拷贝、shell 调用、数学常数、单元转换等）；`qe_base` 类只保留必要的公共接口或被替换为 `StructureAnalyzer`。
- `qebin.py` 提供路径/环境信息，应配合配置文件（如 YAML/JSON）统一维护。
- 编写 `qe_logo`、`qe_logging` 等轻量模块时，整理到 `core` 层，避免循环依赖。

## 测试与验证
- 引入基本单元测试：对参数解析、输入文件生成、lambda 计算等核心功能构建 pytest。
- 为关键工作流准备最小化示例数据 (POSCAR, dyn0, gam.lines 等)，在 CI 中验证模板输出保持一致。
- 编写重构迁移文档：列出旧命令与新命令映射、参数变化、兼容性注意事项。

## 实施步骤建议
1. 建立新目录骨架与 `BaseParameters`/`Workflow` 抽象，迁移 `relax/scf` 作为 POC。
2. 重写输入/提交 writer 提供模板化接口，完成 relax/scf 测试。
3. 按模块分批迁移 `phonon → electron → superconduct → batch/sctk`，逐步删除旧分支。
4. 拆分数据处理逻辑到 `postprocess`，确保绘图/导出脚本统一调用。
5. 调整 CLI 与日志体系，更新文档与示例。
6. 清理遗留函数、移除 `kwargs` + `hasattr` 模式，引入类型检查与测试。

