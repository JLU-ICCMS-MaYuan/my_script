# QE计算工具包

一个模块化、高效的Quantum ESPRESSO (QE)计算工具包，提供完整的工作流管理、批量计算和结果分析功能。

## 主要特性

- ✅ **模块化架构**: 高内聚、低耦合，易于维护和扩展
- ✅ **完整工作流**: 支持Relax、SCF、Phonon、Electron、Superconductivity计算
- ✅ **批量处理**: Pipeline系统支持多结构批量计算
- ✅ **实时监控**: Rich库提供美观的实时进度显示
- ✅ **智能并行**: 混合并行策略（结构级+步骤级）
- ✅ **Publication绘图**: matplotlib + plotly生成发表级图表
- ✅ **工具集成**: 重构整合qe_toolkit常用工具

## 项目结构

```
qe/
├── analysis/           # 数据分析模块
│   ├── plotting/      # 绘图（能带、DOS、声子）
│   ├── superconductivity/  # 超导Tc计算
│   └── thermodynamics/     # 热力学分析
├── config/            # 配置类（Relax/SCF/Phonon/Electron/SC）
├── core/              # 核心模块（常数、异常、类型）
├── io/                # 输入输出
│   ├── readers/      # 输出文件读取器
│   └── writers/      # 输入文件生成器
├── pipelines/         # Pipeline系统
│   ├── base.py       # 基类
│   └── templates/    # 预定义模板
├── scheduler/         # 任务调度
│   ├── task.py       # 任务定义
│   └── executor.py   # 混合并行执行器
├── utils/             # 工具模块
│   ├── checkers/     # 任务状态检查
│   ├── progress.py   # 进度监控
│   └── structure.py  # 结构处理
├── workflows/         # 单步工作流
│   ├── relax.py      # 结构优化
│   ├── scf.py        # 自洽计算
│   ├── phonon.py     # 声子计算
│   ├── electron.py   # 电子性质
│   └── superconductivity.py  # 超导计算
└── examples/          # 使用示例
```

## 快速开始

### 1. 安装依赖

```bash
pip install numpy matplotlib plotly rich
```

### 2. 基本使用

#### 方式一：使用预定义Pipeline模板

```python
from pathlib import Path
from pipelines.templates import RelaxElectronPipeline

# 准备结构文件
structures = list(Path("./structures").glob("*.vasp"))

# 配置参数
config = {
    'ecutwfc': 80,
    'ecutrho': 640,
    'kpoints': (8, 8, 8),
    'pseudopotentials': {
        'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
    },
    'nprocs': 8,
    'qe_bin': '/path/to/qe/bin/',
}

# 创建Pipeline
pipeline = RelaxElectronPipeline(
    structures=structures,
    work_dir=Path("./calculations"),
    config=config,
)

# 运行（带实时监控）
pipeline.run_with_monitor(max_workers=4)
```

#### 方式二：单步工作流

```python
from pathlib import Path
from workflows.scf import ScfWorkflow
from config.base import BaseConfig

# 创建配置
config = BaseConfig()
config.set_param('ecutwfc', 80)
config.set_param('ecutrho', 640)
config.set_param('kpoints', (8, 8, 8))

# 创建工作流
workflow = ScfWorkflow(
    config=config,
    structure_file=Path("structure.vasp")
)

# 运行
workflow.run(submit=False)  # 仅生成输入文件
```

### 3. 使用任务检查器

#### 检查声子计算进度

```bash
# 检查q点1-10的声子计算状态
python qe/utils/checkers/phonon_checker.py 1 10

# 输出漂亮的表格，显示：
# - DYN文件状态（完成/虚频/未完成）
# - ELPH目录状态
# - 不可约表示进度
```

#### 检查任务完成情况

```bash
# 递归检查relax和scf任务
python qe/utils/checkers/relax_checker.py -d 2 -t relax scf

# 显示：
# - 成功/失败/无输出统计
# - 成功率
# - 失败任务列表
```

### 4. 超导Tc计算

```python
from pathlib import Path
from analysis.superconductivity.tc_calculator import TcCalculator

calc = TcCalculator()

# 读取alpha2F数据
alpha2f_file = Path("_ph0/fildyn.a2Fsave")
frequencies, alpha2f = calc.read_alpha2f(alpha2f_file)

# 计算λ和ω_log
lambda_val, omega_log = calc.calculate_lambda(frequencies, alpha2f)

# 计算Tc（多个μ*值）
for mu in [0.10, 0.13, 0.15]:
    tc_mcm = calc.mcmillan_tc(lambda_val, omega_log, mu)
    tc_ad = calc.allen_dynes_tc(lambda_val, omega_log, mu)

    print(f"μ* = {mu:.2f}:")
    print(f"  Tc (McMillan)    = {tc_mcm:6.2f} K")
    print(f"  Tc (Allen-Dynes) = {tc_ad:6.2f} K")
```

## 可用的Pipeline模板

### 1. RelaxPhonoSuperconductivityPipeline
完整的超导性质计算流程：

Relax → SCF → Phonon → Q2R → Matdyn → PhononDOS → El-Ph → Lambda → Tc

### 2. RelaxElectronPipeline
结构优化 + 电子性质：

Relax → SCF → Electron (Bands+DOS)

### 3. ScfPhononPipeline
声子谱计算：

SCF → Phonon → Q2R → Matdyn → PhononDOS

### 4. ScfPhononElphPipeline
声子 + 电声耦合：

SCF → Phonon → El-Ph → Q2R → Matdyn → PhononDOS

## 实时进度监控

所有Pipeline都支持Rich库的实时进度监控：

```python
from utils.progress import ProgressMonitor

# 创建监控器
monitor = ProgressMonitor(total_tasks=len(tasks))

# 运行Pipeline（带监控）
monitor.run_with_monitor(executor.execute, tasks)

# 显示：
# - 实时进度条
# - 任务状态表格（结构/步骤/状态/进度/耗时）
# - 统计信息（成功/失败/运行/等待）
# - ETA（预计剩余时间）
```

## 命令行工具

### 声子任务检查
```bash
python -m qe.utils.checkers.phonon_checker 1 10
```

### 任务状态检查
```bash
python -m qe.utils.checkers.relax_checker -d 2 --show-failed 20
```

### Lambda文件备份
```bash
python -m qe.analysis.superconductivity.lambda_backup 0.10 0.13 0.15
```

## 特性详解

### 1. 混合并行策略

- **结构级并行**: 多个结构同时计算
- **步骤串行**: 每个结构的步骤按顺序执行
- **自动错误处理**: 失败结构跳过，不影响其他结构

### 2. Publication-Quality绘图

- **静态图**: matplotlib，300 DPI PNG/PDF
- **交互图**: Plotly HTML
- **批量绘图**: 自动生成HTML报告

### 3. 模块化设计

- 每个workflow独立，易于测试
- 配置与执行分离
- 清晰的依赖关系

## 配置示例

```python
config = {
    # 基本参数
    'ecutwfc': 80,          # 波函数截断能（Ry）
    'ecutrho': 640,         # 电荷密度截断能（Ry）
    'kpoints': (8, 8, 8),   # K点网格
    'qpoints': '6 6 6',     # Q点网格

    # 收敛参数
    'conv_thr': '1.0d-8',
    'mixing_beta': 0.3,

    # 占据
    'occupations': 'smearing',
    'smearing': 'mp',
    'degauss': 0.02,

    # 赝势
    'pseudopotentials': {
        'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
        'S': 'S.pbe-n-rrkjus_psl.1.0.0.UPF',
    },

    # 并行
    'nprocs': 8,
    'qe_bin': '~/soft/qe-7.4.1/bin/',
}
```

## 测试

```bash
# 运行单元测试
pytest tests/unit/

# 运行集成测试
pytest tests/integration/

# 使用H3S测试结构
python examples/example_with_progress.py
```

## 开发规范

- 使用type hints
- 完整的docstring（Google风格）
- 模块化设计
- 单元测试覆盖

## 许可证

MIT License

## 作者

Claude (2025-11-20)

## 致谢

重构自原有的QE工具包和qe_toolkit工具集。
