# Repository Guidelines
本指南为贡献者提供在`my_script`仓库协作的核心约定，涵盖目录结构、开发流程、风格约束与提交流程，优先保持简洁、可复现与可追踪。

## 项目结构与模块组织
- `qe/`、`vasp/`、`epw/`、`structuregenerator/`：核心CLI模块，分别负责QE、VASP、EPW与结构生成任务，可独立开发与发布。
- `common/`、`mytoolkit/`：通用工具、脚本与配置，新增公共逻辑优先集中维护。
- `example/`、`notebook/`：示例与实验性脚本，便于快速复现计算流程。
- `test/`：测试用例建议放置于此（命名`test_*.py`），与被测模块平行。
- `pyproject.toml`、`requirements.txt`：依赖与工具链声明，变更需同步说明原因。

## 构建、测试与开发命令
- 创建环境：`python -m venv .venv && source .venv/bin/activate`，然后`pip install -r requirements.txt`或`pip install -e .[dev]`。
- 本地运行：安装后可直接执行`qe --help`、`vasp --help`、`epw --help`验证CLI是否正常加载。
- 格式化与检查：`black .`（行宽88）与`ruff check .`（行宽100）保持一致风格。
- 测试：`pytest -v --cov=. --cov-report=term-missing`，提交前确保核心逻辑覆盖≥80%。
- 构建发行物（如需发布）：`python -m build`生成sdist与wheel，产物存放于`dist/`。

## 编码风格与命名规范
- 语言：Python 3.7+，缩进4空格，鼓励完整类型注解与显式错误处理。
- 命名：模块与文件使用snake_case，类用PascalCase，函数/变量用snake_case，常量全大写。
- CLI相关脚本保持与模块同名（如`qe_main.py`），输入输出路径使用绝对或仓库相对路径并提供默认值。
- 提交前运行`black`与`ruff`；需要静态检查时执行`mypy .`保证接口一致性。

## 测试指南
- 框架：Pytest，测试文件命名`test_*.py`，针对CLI入口至少覆盖参数解析、配置优先级与核心流程分支。
- 数据：大型结构或外部势文件使用精简示例或mock，避免提交敏感数据。
- 断言：优先验证关键输出（生成文件、stdout片段、返回码），对计算耗时逻辑拆分可测的纯函数。

## 提交与PR规范
- 提交信息使用Conventional Commits，例如`feat: 支持qe声子批量分割`、`fix: 修正vasp kspacing默认值`，保持单次提交聚焦单一变更。
- PR需描述变更目的、核心修改点、测试结果（命令及摘要），涉及接口变更时补充升级指引；链接相关Issue/讨论。
- 对新增脚本或配置，附简短使用示例或路径说明，便于审核与回溯。
