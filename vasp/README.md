# VASP 计算工具 - 使用手册

一个现代化的VASP第一性原理计算命令行工具，提供完整的计算工作流和批量处理功能。

## 目录

- [快速开始](#快速开始)
- [安装配置](#安装配置)
- [基本使用](#基本使用)
- [详细教程](#详细教程)
- [JSON配置文件](#json配置文件)
- [批量计算](#批量计算)
- [高级功能](#高级功能)
- [常见问题](#常见问题)

---

## 快速开始

### 查看帮助信息

```bash
# 查看主命令帮助
vasp --help

# 查看子命令帮助
vasp relax --help
vasp electronic --help
vasp phonon --help
vasp md --help
```

### 最简单的例子

```bash
# 1. 结构优化
vasp relax -i POSCAR -w ./relax --kspacing 0.2 --potcar-dir ~/pot/vasp_pot/potpaw_PBE54

# 2. 电子性质计算（完整流程）
vasp electronic -i POSCAR -w ./electronic --kspacing 0.2 --potcar-dir ~/pot/vasp_pot/potpaw_PBE54

# 3. 声子计算
vasp phonon -i POSCAR -w ./phonon --supercell 2 2 2 --kspacing 0.2 --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

---

## 安装配置

### 1. 安装my_script包

```bash
cd /path/to/my_script
pip install -e .
```

安装成功后，`vasp` 命令会自动添加到系统PATH。

### 2. 准备POTCAR库

确保你有VASP的POTCAR赝势库，例如：

```bash
~/pot/vasp_pot/
├── potpaw_PBE54/
│   ├── H/
│   ├── Li/
│   ├── C/
│   └── ...
├── potpaw_LDA/
└── potpaw_GGA/
```

### 3. 配置VASP执行命令

确保VASP可执行文件已正确配置。如果需要MPI并行：

```bash
# 加载Intel OneAPI环境（如果使用）
source ~/intel/oneapi/setvars.sh

# 测试VASP是否可用
mpirun -np 4 vasp_std
```

---

## 基本使用

### 1. 结构优化 (relax)

#### 单个结构优化

```bash
vasp relax \
  -i POSCAR \
  -w ./calc/relax \
  --kspacing 0.2 \
  --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54 \
  --potcar-type PBE \
  -j bash
```

**参数说明：**
- `-i, --input`: 输入结构文件（POSCAR格式）
- `-w, --work-dir`: 工作目录
- `--kspacing`: K点间距（埃的倒数）
- `--encut`: 截断能（eV）
- `--potcar-dir`: POTCAR库路径
- `--potcar-type`: POTCAR类型（PBE/LDA/PW91）
- `-j, --job-system`: 队列系统（bash/slurm/pbs）

#### 批量结构优化

```bash
# 方法1：使用--batch标志
vasp relax \
  -i ./structures/ \
  -w ./batch_relax \
  --batch \
  --kspacing 0.2 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54

# 方法2：自动检测（目录中包含多个.vasp或POSCAR文件）
vasp relax \
  -i ./structures/ \
  -w ./batch_relax \
  --kspacing 0.2 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

### 2. 电子性质计算 (electronic)

完整的电子性质计算流程包括：
1. **Relax**: 结构优化
2. **SCF**: 自洽计算
3. **DOS**: 电子态密度
4. **Band**: 能带结构
5. **ELF** (可选): 电子局域函数
6. **COHP** (可选): 晶体轨道哈密顿布居
7. **Plotting**: 自动绘图

#### 基本用法

```bash
vasp electronic \
  -i POSCAR \
  -w ./electronic \
  --kspacing 0.2 \
  --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

#### 包含ELF和COHP分析

```bash
vasp electronic \
  -i POSCAR \
  -w ./electronic \
  --kspacing 0.2 \
  --encut 500 \
  --include-elf \
  --include-cohp \
  --dos-type element_spd \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

**参数说明：**
- `--include-elf`: 包含电子局域函数计算
- `--include-cohp`: 包含COHP分析
- `--dos-type`: DOS投影类型（element/spd/element_spd）

### 3. 声子性质计算 (phonon)

完整的声子性质计算流程包括：
1. **Relax**: 结构优化
2. **Phonon**: 声子计算（有限位移法或DFPT）
3. **Band**: 声子谱
4. **DOS**: 声子态密度
5. **Plotting**: 自动绘图

#### 有限位移法

```bash
vasp phonon \
  -i POSCAR \
  -w ./phonon \
  --supercell 2 2 2 \
  --method disp \
  --kspacing 0.2 \
  --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

#### DFPT方法

```bash
vasp phonon \
  -i POSCAR \
  -w ./phonon \
  --supercell 2 2 2 \
  --method dfpt \
  --kspacing 0.2 \
  --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

**参数说明：**
- `--supercell X Y Z`: 超胞大小
- `--method`: 计算方法（disp=有限位移法, dfpt=密度泛函微扰理论）

### 4. 分子动力学 (md)

```bash
vasp md \
  -i POSCAR \
  -w ./md \
  --kspacing 0.3 \
  --temp 300 \
  --time 1000 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

---

## 详细教程

### 教程1：LiH电子性质完整计算

这是一个从头到尾的完整示例。

#### 第1步：准备结构文件

创建 `LiH.vasp` 文件：

```
Li1H1
1.0
  2.84  0.00  0.00
  0.00  2.84  0.00
  0.00  0.00  2.84
Li H
1 1
direct
  0.0  0.0  0.0  Li
  0.5  0.5  0.5  H
```

#### 第2步：创建JSON配置文件

创建 `vasp_config.json`：

```json
{
  "potcar_dir": "~/pot/vasp_pot/potpaw_PBE54",
  "potcar_type": "PBE",
  "kspacing": 0.2,
  "encut": 500,
  "job_system": "bash",
  "include_elf": false,
  "include_cohp": false,
  "dos_type": "element"
}
```

#### 第3步：运行电子性质计算

```bash
vasp electronic \
  -i LiH.vasp \
  -w ./LiH_electronic \
  --json vasp_config.json
```

#### 第4步：查看结果

计算完成后，目录结构：

```
LiH_electronic/
├── 01_relax/
│   ├── POSCAR
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   ├── CONTCAR
│   └── OUTCAR
├── 02_scf/
│   └── ...
├── 03_dos/
│   └── ...
├── 04_band/
│   └── ...
└── plots/
    ├── dos.png
    └── band.png
```

### 教程2：MgB2声子谱计算

超导材料MgB2的声子谱计算示例。

#### 第1步：准备MgB2结构

创建 `MgB2.vasp`：

```
MgB2
1.0
  3.08  0.00  0.00
 -1.54  2.67  0.00
  0.00  0.00  3.52
Mg B
1 2
direct
  0.0000  0.0000  0.0000  Mg
  0.3333  0.6667  0.5000  B
  0.6667  0.3333  0.5000  B
```

#### 第2步：声子计算

```bash
vasp phonon \
  -i MgB2.vasp \
  -w ./MgB2_phonon \
  --supercell 2 2 2 \
  --method disp \
  --kspacing 0.15 \
  --encut 600 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

#### 第3步：分析虚频

声子计算完成后，检查是否有虚频：

```bash
cd MgB2_phonon/phonon
grep "f  =" phonon_band.dat | awk '{if($3<0) print}'
```

### 教程3：批量结构优化

假设你有一批氢化物结构需要优化。

#### 第1步：准备结构目录

```
structures/
├── LiH.vasp
├── NaH.vasp
├── KH.vasp
├── RbH.vasp
└── CsH.vasp
```

#### 第2步：批量优化

```bash
vasp relax \
  -i ./structures/ \
  -w ./batch_relax \
  --batch \
  --parallel \
  --max-workers 5 \
  --kspacing 0.2 \
  --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

**参数说明：**
- `--parallel`: 启用并行批量计算
- `--max-workers 5`: 最多同时运行5个任务

#### 第3步：检查结果

```bash
cd batch_relax
ls -d */

# 输出：
# LiH/
# NaH/
# KH/
# RbH/
# CsH/

# 查看每个计算的能量
for dir in */; do
  echo "$dir:"
  grep "energy  without" $dir/01_relax/OUTCAR | tail -1
done
```

### 教程4：高精度电子结构计算

对于需要发表的高精度计算。

#### 第1步：创建高精度配置

创建 `high_precision.json`：

```json
{
  "potcar_dir": "~/pot/vasp_pot/potpaw_PBE54",
  "potcar_type": "PBE",
  "kspacing": 0.1,
  "encut": 800,
  "job_system": "slurm",
  "include_elf": true,
  "include_cohp": true,
  "dos_type": "element_spd"
}
```

#### 第2步：运行计算

```bash
vasp electronic \
  -i structure.vasp \
  -w ./high_precision \
  --json high_precision.json \
  --log-level DEBUG
```

---

## JSON配置文件

### 为什么使用JSON配置？

JSON配置文件可以：
1. 存储常用参数，避免重复输入
2. 团队共享配置，保证计算一致性
3. 记录计算参数，便于后期复现

### 配置优先级

```
命令行参数 > JSON配置 > 代码默认值
```

### 完整配置示例

```json
{
  "_comment": "VASP计算通用配置",

  "potcar_dir": "~/pot/vasp_pot/potpaw_PBE54",
  "potcar_type": "PBE",

  "kspacing": 0.2,
  "encut": 500,

  "job_system": "bash",
  "max_workers": 4,
  "parallel": false,

  "include_elf": false,
  "include_cohp": false,
  "dos_type": "element",

  "supercell": [2, 2, 2],
  "method": "disp",

  "log_level": "INFO"
}
```

### 使用配置文件

```bash
# 使用默认配置
vasp electronic -i POSCAR -w ./calc --json config.json

# 覆盖部分配置
vasp electronic -i POSCAR -w ./calc --json config.json --kspacing 0.15 --encut 600
```

---

## 批量计算

### 自动批量检测

VASP CLI会自动检测是否需要批量模式：

```bash
# 输入是目录且包含多个结构文件 → 自动启用批量模式
vasp relax -i ./structures/ -w ./batch
```

检测规则：
- 目录中存在 `*.vasp` 文件
- 目录中存在 `*.POSCAR` 文件
- 目录中存在 `POSCAR*` 文件

### 并行批量计算

```bash
# 并行处理4个结构
vasp electronic \
  -i ./structures/ \
  -w ./batch \
  --parallel \
  --max-workers 4 \
  --kspacing 0.2 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54
```

**注意：**
- 每个worker会占用1个计算任务
- 确保有足够的计算资源
- 建议 `max-workers` ≤ CPU核心数

### 批量计算目录结构

```
batch/
├── structure1/
│   ├── 01_relax/
│   ├── 02_scf/
│   └── ...
├── structure2/
│   ├── 01_relax/
│   ├── 02_scf/
│   └── ...
└── batch_summary.json
```

---

## 高级功能

### 1. 自定义VASP执行命令

如果需要特殊的MPI设置：

```bash
# 使用Intel MPI
source ~/intel/oneapi/setvars.sh
vasp electronic -i POSCAR -w ./calc --json config.json

# 或在JSON中指定
{
  "execmd": "mpirun -np 8 vasp_std",
  ...
}
```

### 2. 调试模式

```bash
vasp electronic \
  -i POSCAR \
  -w ./calc \
  --json config.json \
  --log-level DEBUG
```

这会输出详细的调试信息。

### 3. 仅生成输入文件

目前CLI主要用于准备和提交任务。如果只想生成输入文件而不运行：

```python
# 使用Python API
from vasp.pipelines import ElectronicPropertiesPipeline

pipeline = ElectronicPropertiesPipeline(
    structure_file="POSCAR",
    work_dir="./calc",
    kspacing=0.2,
    encut=500
)

# 只准备输入文件
pipeline.prepare_inputs()
```

### 4. 指定DOS投影类型

```bash
# 投影到元素
vasp electronic -i POSCAR -w ./calc --dos-type element

# 投影到spd轨道
vasp electronic -i POSCAR -w ./calc --dos-type spd

# 同时投影元素和spd
vasp electronic -i POSCAR -w ./calc --dos-type element_spd
```

---

## 常见问题

### Q1: 如何检查计算是否完成？

```bash
# 检查OUTCAR
grep "reached required accuracy" OUTCAR

# 检查是否正常结束
tail OUTCAR
```

### Q2: POTCAR找不到怎么办？

确保：
1. `--potcar-dir` 路径正确
2. 目录中存在对应元素的POTCAR
3. 使用正确的 `--potcar-type`（PBE/LDA/PW91）

示例目录结构：
```
~/pot/vasp_pot/potpaw_PBE54/
├── H/
│   └── POTCAR
├── Li/
│   └── POTCAR
└── ...
```

### Q3: 如何设置K点间距？

K点间距与精度的关系：
- `kspacing = 0.5`: 粗糙，快速测试
- `kspacing = 0.3`: 中等精度
- `kspacing = 0.2`: 常规收敛精度
- `kspacing = 0.15`: 高精度
- `kspacing = 0.1`: 发表级精度

### Q4: 批量计算失败了怎么办？

批量计算中单个结构失败不会影响其他结构：

```bash
# 查看失败的任务
cat batch/batch_summary.json

# 单独重新运行失败的结构
vasp electronic -i failed_structure.vasp -w ./retry --json config.json
```

### Q5: 如何选择截断能？

一般规则：
```bash
# 查看POTCAR推荐的ENMAX
grep ENMAX ~/pot/vasp_pot/potpaw_PBE54/*/POTCAR

# 使用 1.3 × ENMAX_max
vasp relax -i POSCAR -w ./calc --encut 500
```

### Q6: 声子谱出现虚频怎么办？

1. 检查结构是否充分优化
2. 增大超胞尺寸
3. 提高K点密度
4. 检查是否是力学不稳定

```bash
# 重新优化结构
vasp relax -i POSCAR -w ./rerelax --kspacing 0.15 --encut 600

# 使用更大超胞
vasp phonon -i CONTCAR -w ./phonon_large --supercell 3 3 3
```

### Q7: 内存不足怎么办？

减少内存占用：
1. 减小K点密度
2. 减少截断能
3. 使用更多节点并行

### Q8: 如何在超算上运行？

```bash
# 使用Slurm队列系统
vasp electronic \
  -i POSCAR \
  -w ./calc \
  --json config.json \
  -j slurm

# 使用PBS队列系统
vasp electronic \
  -i POSCAR \
  -w ./calc \
  --json config.json \
  -j pbs
```

---

## 架构说明

VASP模块采用模块化架构，详细信息参见[架构文档](./ARCHITECTURE.md)。

主要组件：
- `workflows/`: 单步工作流（relax, phonon, electron, md）
- `pipelines/`: 多步骤自动化流程
- `io/`: 输入输出文件处理
- `analysis/`: 数据分析和绘图
- `cli.py`: 命令行接口

---

## 版本信息

- **当前版本**: 2.0.0
- **CLI版本**: 1.0.0
- **更新日期**: 2025-11-20
- **开发工具**: Claude

---

## 许可证

MIT License

---

## 相关资源

- [VASP官方文档](https://www.vasp.at/wiki/index.php/The_VASP_Manual)
- [VASP论坛](https://www.vasp.at/forum/)
- [Materials Project](https://materialsproject.org/)

---

## 技术支持

如遇问题，请：
1. 查看本文档的[常见问题](#常见问题)章节
2. 检查 `--log-level DEBUG` 输出
3. 在项目Issues中提问

---

**祝你计算顺利！** 🚀
