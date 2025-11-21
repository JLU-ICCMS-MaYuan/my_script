# Quantum ESPRESSO (QE) 计算工具 - 使用手册

一个现代化的Quantum ESPRESSO第一性原理计算命令行工具，提供完整的计算工作流和超导性质分析功能。

## 目录

- [快速开始](#快速开始)
- [安装配置](#安装配置)
- [基本使用](#基本使用)
- [详细教程](#详细教程)
- [JSON配置文件](#json配置文件)
- [超导计算专题](#超导计算专题)
- [声子计算专题](#声子计算专题)
- [常见问题](#常见问题)

---

## 快速开始

### 查看帮助信息

```bash
# 查看主命令帮助
qe --help

# 查看子命令帮助
qe relax --help
qe scf --help
qe phonon --help
qe electron --help
qe superconductivity --help
```

### 最简单的例子

```bash
# 1. 结构优化
qe relax -i structure.cif -w ./relax --ecutwfc 80 --kpoints 16 16 16 --pp-dir ~/pot/qe_pp

# 2. 自洽计算
qe scf -i relax.out -w ./scf --ecutwfc 80 --kpoints 16 16 16 --pp-dir ~/pot/qe_pp

# 3. 声子计算
qe phonon -i scf.out -w ./phonon --qpoints 4 4 4 --pp-dir ~/pot/qe_pp

# 4. 超导性质计算
qe superconductivity -i relax.out -w ./sc --method McAD --mu-star 0.1 --pp-dir ~/pot/qe_pp
```

---

## 安装配置

### 1. 安装my_script包

```bash
cd /path/to/my_script
pip install -e .
```

安装成功后，`qe` 命令会自动添加到系统PATH。

### 2. 准备QE赝势库

确保你有QE的赝势文件，例如：

```bash
~/pot/qe_pp/
├── H.pbe-rrkjus_psl.1.0.0.UPF
├── Li.pbe-sl-rrkjus_psl.1.0.0.UPF
├── C.pbe-n-rrkjus_psl.1.0.0.UPF
├── S.pbe-n-rrkjus_psl.1.0.0.UPF
└── ...
```

推荐使用PSLibrary或SSSP赝势库：
- **PSLibrary**: http://pseudopotentials.quantum-espresso.org/
- **SSSP**: https://www.materialscloud.org/discover/sssp

### 3. 配置QE执行命令

确保QE可执行文件已正确配置：

```bash
# 测试pw.x是否可用
which pw.x
mpirun -np 4 pw.x --version

# 测试ph.x是否可用
which ph.x
```

---

## 基本使用

### 1. 结构优化 (relax)

#### 单个结构优化

```bash
qe relax \
  -i structure.cif \
  -w ./relax \
  --ecutwfc 80 \
  --ecutrho 640 \
  --kpoints 16 16 16 \
  --pp-dir ~/pot/qe_pp \
  --mode vc-relax \
  -j bash
```

**参数说明：**
- `-i, --input`: 输入结构文件（支持cif/vasp/xsf等格式）
- `-w, --work-dir`: 工作目录
- `--ecutwfc`: 波函数截断能（Ry）
- `--ecutrho`: 电荷密度截断能（Ry，通常 = 8 × ecutwfc）
- `--kpoints X Y Z`: K点网格
- `--pp-dir`: 赝势文件目录
- `--mode`: 优化模式（relax固定晶胞/vc-relax变胞）
- `-j, --job-system`: 队列系统（bash/slurm/pbs）

#### 变胞优化（推荐用于结构搜索）

```bash
qe relax \
  -i structure.cif \
  -w ./relax \
  --mode vc-relax \
  --press 0.0 \
  --ecutwfc 80 \
  --kpoints 16 16 16 \
  --pp-dir ~/pot/qe_pp
```

### 2. 自洽计算 (scf)

```bash
qe scf \
  -i relax.out \
  -w ./scf \
  --ecutwfc 80 \
  --ecutrho 640 \
  --kpoints 16 16 16 \
  --degauss 0.02 \
  --pp-dir ~/pot/qe_pp
```

**重要参数：**
- `--degauss`: 展宽（Ry），通常0.01-0.05
  - 金属体系：建议0.02
  - 半导体/绝缘体：可以更小
- `--occupations`: 占据方式
  - `smearing`: 金属（默认）
  - `fixed`: 绝缘体

### 3. 声子计算 (phonon)

声子计算是QE的核心功能，支持多种模式。

#### 不分q点计算（小体系）

```bash
qe phonon \
  -i scf.out \
  -w ./phonon \
  --qpoints 4 4 4 \
  --tr2-ph 1.0e-14 \
  --pp-dir ~/pot/qe_pp
```

#### 分割q点计算（大体系推荐）

```bash
# 第1步：生成q点分割列表
qe phonon \
  -i scf.out \
  -w ./phonon \
  --qpoints 4 4 4 \
  --split \
  --pp-dir ~/pot/qe_pp

# 第2步：提交各个q点计算（手动或使用队列）
# 每个q点目录会生成独立的计算任务

# 第3步：合并q点结果
qe phonon \
  -i scf.out \
  -w ./phonon \
  --merge \
  --pp-dir ~/pot/qe_pp
```

#### 计算电声耦合常数

```bash
qe phonon \
  -i scf.out \
  -w ./phonon \
  --qpoints 6 6 6 \
  --compute-lambda \
  --pp-dir ~/pot/qe_pp
```

**参数说明：**
- `--qpoints X Y Z`: q点网格
- `--tr2-ph`: 声子收敛阈值
- `--split`: 分割q点计算
- `--merge`: 合并q点结果
- `--compute-lambda`: 计算λ和ωlog

### 4. 电子性质计算 (electron)

#### 能带计算

```bash
qe electron \
  -i scf.out \
  -w ./band \
  --mode band \
  --nbnd 50 \
  --pp-dir ~/pot/qe_pp
```

#### 态密度计算

```bash
qe electron \
  -i scf.out \
  -w ./dos \
  --mode dos \
  --kpoints 24 24 24 \
  --pp-dir ~/pot/qe_pp
```

#### 能带+DOS联合计算

```bash
qe electron \
  -i scf.out \
  -w ./electronic \
  --mode band+dos \
  --nbnd 50 \
  --kpoints 24 24 24 \
  --pp-dir ~/pot/qe_pp
```

### 5. 超导性质计算 (superconductivity)

#### McMillan-Allen-Dynes方法

```bash
qe superconductivity \
  -i scf.out \
  -w ./sc_mcad \
  --method McAD \
  --mu-star 0.1 0.13 0.15 \
  --pp-dir ~/pot/qe_pp
```

#### Eliashberg方程方法

```bash
qe superconductivity \
  -i scf.out \
  -w ./sc_eliashberg \
  --method eliashberg \
  --mu-star 0.1 \
  --ntemp 50 \
  --pp-dir ~/pot/qe_pp
```

#### q2r方法（从动力学矩阵）

```bash
qe superconductivity \
  -i scf.out \
  -w ./sc_q2r \
  --method q2r \
  --qpoints 6 6 6 \
  --pp-dir ~/pot/qe_pp
```

**参数说明：**
- `--method`: 计算方法
  - `McAD`: McMillan-Allen-Dynes公式（快速估算）
  - `eliashberg`: 各向异性Eliashberg方程（精确）
  - `q2r`: 从动力学矩阵计算
- `--mu-star`: 库仑赝势μ*（可以指定多个值）
- `--ntemp`: 温度点数（用于Eliashberg方程）

---

## 详细教程

### 教程1：H3S超导性质完整计算

H3S是著名的高温超导材料。这个教程展示完整的超导计算流程。

#### 第1步：准备H3S结构

创建 `H3S.cif`：

```
data_H3S
_cell_length_a    3.1
_cell_length_b    3.1
_cell_length_c    3.1
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'Im-3m'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
S  0.0  0.0  0.0
H  0.5  0.5  0.0
```

#### 第2步：创建QE配置文件

创建 `qe_config.json`：

```json
{
  "pp_dir": "~/pot/qe_pp",
  "ecutwfc": 80,
  "ecutrho": 640,
  "kpoints": [16, 16, 16],
  "qpoints": [6, 6, 6],
  "degauss": 0.02,
  "occupations": "smearing",
  "job_system": "bash",
  "mu_star": [0.10, 0.13, 0.15]
}
```

#### 第3步：结构优化

```bash
qe relax \
  -i H3S.cif \
  -w ./01_relax \
  --json qe_config.json \
  --mode vc-relax \
  --press 150  # 150 GPa压强
```

#### 第4步：自洽计算

```bash
qe scf \
  -i 01_relax/relax.out \
  -w ./02_scf \
  --json qe_config.json
```

#### 第5步：声子+电声耦合计算

```bash
qe phonon \
  -i 02_scf/scf.out \
  -w ./03_phonon \
  --json qe_config.json \
  --compute-lambda
```

#### 第6步：计算超导温度

```bash
qe superconductivity \
  -i 02_scf/scf.out \
  -w ./04_sc \
  --method McAD \
  --json qe_config.json
```

#### 第7步：分析结果

```bash
# 查看λ和ωlog
cat 03_phonon/lambda.out

# 查看Tc
cat 04_sc/tc_results.dat
```

预期结果（在150 GPa下）：
- λ ≈ 2.0-2.5
- ωlog ≈ 1500-2000 K
- Tc ≈ 200 K (μ* = 0.1)

### 教程2：MgB2声子谱计算

MgB2是经典的两能隙超导体。

#### 第1步：准备MgB2结构

```bash
# MgB2.cif
cat > MgB2.cif << EOF
data_MgB2
_cell_length_a    3.086
_cell_length_b    3.086
_cell_length_c    3.524
_cell_angle_alpha 90
_cell_angle_beta  90
_cell_angle_gamma 120
_symmetry_space_group_name_H-M 'P6/mmm'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Mg 0.0    0.0    0.0
B  0.3333 0.6667 0.5
EOF
```

#### 第2步：计算流程

```bash
# 优化
qe relax -i MgB2.cif -w ./relax --ecutwfc 80 --kpoints 24 24 18 --pp-dir ~/pot/qe_pp

# 自洽
qe scf -i relax/relax.out -w ./scf --ecutwfc 80 --kpoints 24 24 18 --pp-dir ~/pot/qe_pp

# 声子（分q点计算）
qe phonon -i scf/scf.out -w ./phonon --qpoints 6 6 6 --split --pp-dir ~/pot/qe_pp
```

#### 第3步：检查声子谱

```bash
# 合并结果
qe phonon -i scf/scf.out -w ./phonon --merge --pp-dir ~/pot/qe_pp

# 查看声子DOS
cat phonon/phonon.dos

# 检查E2g模式（约600 cm-1）
grep "E2g" phonon/phonon_modes.dat
```

### 教程3：Li金属Fermi面计算

计算金属Fermi面和态密度。

#### 第1步：准备Li BCC结构

```python
# 用ASE生成
from ase.build import bulk
from ase.io import write

li = bulk('Li', 'bcc', a=3.49)
write('Li.cif', li)
```

#### 第2步：高密度K点SCF

```bash
qe scf \
  -i Li.cif \
  -w ./scf \
  --ecutwfc 60 \
  --kpoints 32 32 32 \
  --degauss 0.02 \
  --occupations smearing \
  --pp-dir ~/pot/qe_pp
```

#### 第3步：计算DOS

```bash
qe electron \
  -i scf/scf.out \
  -w ./dos \
  --mode dos \
  --kpoints 48 48 48 \
  --pp-dir ~/pot/qe_pp
```

#### 第4步：分析Fermi面

```bash
# 查看费米能级
grep "Fermi" scf/scf.out

# 查看费米面处态密度
grep "EFermi" dos/dos.out
```

### 教程4：批量超导材料筛选

假设你有一批候选超导材料需要快速筛选。

#### 第1步：准备结构目录

```
structures/
├── H3S.cif
├── LaH10.cif
├── YH9.cif
└── CaH6.cif
```

#### 第2步：创建批量计算脚本

```bash
#!/bin/bash
# batch_sc_screen.sh

for cif in structures/*.cif; do
    name=$(basename $cif .cif)
    echo "Processing $name..."

    # 优化
    qe relax -i $cif -w ./$name/relax \
        --ecutwfc 80 --kpoints 16 16 16 \
        --pp-dir ~/pot/qe_pp --json qe_config.json

    # 自洽
    qe scf -i ./$name/relax/relax.out -w ./$name/scf \
        --json qe_config.json

    # 声子+λ
    qe phonon -i ./$name/scf/scf.out -w ./$name/phonon \
        --qpoints 4 4 4 --compute-lambda \
        --json qe_config.json

    # Tc估算
    qe superconductivity -i ./$name/scf/scf.out -w ./$name/sc \
        --method McAD --mu-star 0.1 \
        --json qe_config.json
done
```

#### 第3步：提取结果

```bash
# 提取所有材料的Tc
for dir in */sc; do
    name=$(dirname $dir)
    tc=$(grep "Tc (McMillan)" $dir/tc_results.dat | awk '{print $4}')
    echo "$name: Tc = $tc K"
done > tc_summary.txt

# 排序
sort -t':' -k2 -n -r tc_summary.txt
```

---

## JSON配置文件

### 完整QE配置示例

```json
{
  "_comment": "QE计算通用配置",

  "pp_dir": "~/pot/qe_pp",

  "ecutwfc": 80,
  "ecutrho": 640,

  "kpoints": [16, 16, 16],
  "qpoints": [6, 6, 6],

  "degauss": 0.02,
  "occupations": "smearing",
  "smearing": "mp",

  "conv_thr": "1.0e-8",
  "mixing_beta": 0.3,

  "job_system": "bash",

  "mu_star": [0.10, 0.13, 0.15],
  "ntemp": 50,

  "log_level": "INFO"
}
```

### 赝势配置

如果使用不同元素的特定赝势：

```json
{
  "pp_dir": "~/pot/qe_pp",
  "pseudopotentials": {
    "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
    "S": "S.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Li": "Li.pbe-sl-rrkjus_psl.1.0.0.UPF"
  }
}
```

---

## 超导计算专题

### 理论背景

#### McMillan公式

$$
T_c = \frac{\omega_{log}}{1.2} \exp\left[-\frac{1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)}\right]
$$

#### Allen-Dynes公式（更精确）

$$
T_c = \frac{\omega_{log}}{1.2} \exp\left[-\frac{1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)} \cdot f_1 \cdot f_2\right]
$$

其中：
- λ: 电声耦合常数
- ωlog: 对数平均声子频率
- μ*: 库仑赝势（通常0.10-0.15）

### λ和ωlog的计算

#### 从lambda.out提取

```bash
# 查看lambda.out
cat phonon/lambda.out

# 提取λ
lambda=$(grep "lambda" phonon/lambda.out | awk '{print $3}')

# 提取ωlog
omega_log=$(grep "omega_log" phonon/lambda.out | awk '{print $3}')

echo "λ = $lambda"
echo "ω_log = $omega_log K"
```

#### 从a2F.dos文件计算

```bash
# a2F.dos包含Eliashberg谱函数
# 第1列：频率（meV）
# 第2列：α²F(ω)

# 计算λ = 2∫[α²F(ω)/ω]dω
python calc_lambda.py a2F.dos
```

### 超导参数优化

#### 选择合适的μ*

```bash
# 计算多个μ*值
qe superconductivity \
  -i scf.out \
  -w ./sc \
  --method McAD \
  --mu-star 0.08 0.10 0.12 0.13 0.15 \
  --pp-dir ~/pot/qe_pp

# 查看Tc随μ*的变化
cat sc/tc_vs_mustar.dat
```

通常：
- 简单金属：μ* = 0.10
- 过渡金属：μ* = 0.12-0.13
- 强关联体系：μ* = 0.15

#### q点网格收敛性测试

```bash
for nq in 4 6 8 10; do
    qe phonon \
      -i scf.out \
      -w ./phonon_q${nq} \
      --qpoints $nq $nq $nq \
      --compute-lambda \
      --pp-dir ~/pot/qe_pp

    lambda=$(grep "lambda" ./phonon_q${nq}/lambda.out | awk '{print $3}')
    echo "nq=$nq: λ=$lambda"
done
```

### Eliashberg方程详细计算

对于精确的Tc计算，使用各向异性Eliashberg方程：

```bash
# 第1步：高密度k点SCF
qe scf \
  -i relax.out \
  -w ./scf_dense \
  --kpoints 24 24 24 \
  --ecutwfc 80 \
  --pp-dir ~/pot/qe_pp

# 第2步：声子+电声矩阵元
qe phonon \
  -i scf_dense/scf.out \
  -w ./phonon_elph \
  --qpoints 8 8 8 \
  --compute-lambda \
  --save-dvscf \
  --pp-dir ~/pot/qe_pp

# 第3步：求解Eliashberg方程
qe superconductivity \
  -i scf_dense/scf.out \
  -w ./eliashberg \
  --method eliashberg \
  --mu-star 0.1 \
  --ntemp 100 \
  --temp-range 10 300 \
  --pp-dir ~/pot/qe_pp
```

结果文件：
- `eliashberg_gap_T.dat`: 能隙随温度变化
- `eliashberg_tc.dat`: 临界温度

---

## 声子计算专题

### 虚频处理

#### 检查虚频

```bash
# 从dyn文件检查
grep "freq" phonon/*.dyn* | awk '{if($8<0) print}'

# 从声子DOS检查
head -20 phonon/phonon.dos
```

#### 虚频原因

1. **结构未充分优化**
   ```bash
   # 重新优化，更严格的收敛标准
   qe relax -i POSCAR -w ./relax \
       --forc-conv-thr 1.0e-5 \
       --etot-conv-thr 1.0e-8 \
       --ecutwfc 100
   ```

2. **K点不足**
   ```bash
   # 增加K点密度
   qe scf -i relax.out -w ./scf --kpoints 24 24 24
   ```

3. **动力学不稳定**（真实的虚频）
   ```bash
   # 检查是否为Gamma点虚频
   grep "Gamma" phonon/phonon_modes.dat
   ```

### 分q点并行计算

对于大体系，分q点计算可以显著加速。

#### 自动分割

```bash
# 生成q点列表和计算脚本
qe phonon \
  -i scf.out \
  -w ./phonon \
  --qpoints 6 6 6 \
  --split \
  --pp-dir ~/pot/qe_pp

# 这会创建：
# phonon/
# ├── q_001/
# │   └── ph.in
# ├── q_002/
# │   └── ph.in
# └── ...
```

#### 提交并行任务

```bash
# 使用Slurm数组任务
cat > submit_phonon.sh << 'EOF'
#!/bin/bash
#SBATCH --array=1-216  # 6x6x6 = 216个q点
#SBATCH --ntasks=8

cd phonon/q_$(printf "%03d" $SLURM_ARRAY_TASK_ID)
mpirun -np 8 ph.x -npool 4 < ph.in > ph.out
EOF

sbatch submit_phonon.sh
```

#### 合并结果

```bash
# 等所有q点完成后
qe phonon \
  -i scf.out \
  -w ./phonon \
  --merge \
  --pp-dir ~/pot/qe_pp

# 检查合并结果
ls phonon/*.dyn*  # 应该有完整的dyn文件
```

### 声子DOS和声子谱

#### 计算声子DOS

```bash
# 从动力学矩阵计算DOS
qe phonon \
  -i scf.out \
  -w ./phonon \
  --mode phdos \
  --qpoints-dos 20 20 20 \
  --pp-dir ~/pot/qe_pp
```

#### 绘制声子谱

```bash
# 沿高对称路径
qe phonon \
  -i scf.out \
  -w ./phonon \
  --mode phband \
  --band-path "G-X-M-G-R" \
  --pp-dir ~/pot/qe_pp
```

---

## 常见问题

### Q1: 如何选择截断能？

```bash
# 方法1：查看赝势推荐值
grep "Suggested minimum cutoff" ~/pot/qe_pp/H.pbe*.UPF

# 方法2：收敛性测试
for ecut in 60 70 80 90 100; do
    qe scf -i POSCAR -w ./test_ecut${ecut} \
        --ecutwfc $ecut --ecutrho $((ecut*8)) \
        --kpoints 16 16 16 --pp-dir ~/pot/qe_pp

    energy=$(grep "!" ./test_ecut${ecut}/scf.out | tail -1 | awk '{print $5}')
    echo "$ecut $energy"
done
```

通常：
- 轻元素（H, C, O）：60-80 Ry
- 过渡金属：80-100 Ry
- 含f电子元素：100-120 Ry

### Q2: SCF不收敛怎么办？

```bash
# 1. 减小mixing_beta
qe scf -i POSCAR -w ./scf --mixing-beta 0.1

# 2. 更改混合模式
qe scf -i POSCAR -w ./scf --mixing-mode "plain"

# 3. 增加最大迭代步数
qe scf -i POSCAR -w ./scf --electron-maxstep 200

# 4. 使用更好的初始猜测
cp previous_calc/charge-density.dat ./
qe scf -i POSCAR -w ./scf --startingpot file
```

### Q3: 声子计算中断了怎么续算？

```bash
# QE声子计算支持自动续算
# 只需要重新运行相同的命令
qe phonon -i scf.out -w ./phonon --qpoints 6 6 6 --pp-dir ~/pot/qe_pp

# QE会检测到 _ph0/ 目录并从中断处继续
```

### Q4: 如何检查计算是否完成？

```bash
# 检查SCF
grep "convergence has been achieved" scf/scf.out

# 检查声子
grep "End of self-consistent calculation" phonon/ph.out

# 检查是否正常结束
tail -20 scf/scf.out
```

### Q5: degauss如何设置？

degauss（展宽）的选择很重要：

```bash
# 金属体系
qe scf -i POSCAR -w ./scf --degauss 0.02 --smearing mp

# 半导体/绝缘体
qe scf -i POSCAR -w ./scf --degauss 0.01 --occupations fixed

# 窄带隙材料
qe scf -i POSCAR -w ./scf --degauss 0.005 --smearing gauss
```

### Q6: 内存不足怎么办？

```bash
# 减少对角化缓存
qe scf -i POSCAR -w ./scf --diago-david-ndim 2

# 使用更多的并行池
mpirun -np 32 pw.x -npool 8 < scf.in > scf.out

# 减少K点（如果可行）
qe scf -i POSCAR -w ./scf --kpoints 12 12 12
```

### Q7: 如何设置压强？

```bash
# 在指定压强下优化结构
qe relax \
  -i POSCAR \
  -w ./relax_150GPa \
  --mode vc-relax \
  --press 150 \  # 单位：GPa
  --pp-dir ~/pot/qe_pp
```

### Q8: 如何并行计算？

```bash
# K点池并行（推荐）
mpirun -np 32 pw.x -npool 8 < scf.in > scf.out
# 32个核，分成8个池，每池计算4个k点

# 平面波并行
mpirun -np 32 pw.x -ndiag 4 < scf.in > scf.out

# 图像并行（NEB, PHonon等）
mpirun -np 64 ph.x -npool 16 -nimage 4 < ph.in > ph.out
```

---

## 架构说明

QE模块采用模块化架构：

```
qe/
├── cli.py              # 命令行接口
├── workflows/          # 单步工作流
│   ├── relax.py
│   ├── scf.py
│   ├── phonon.py
│   ├── electron.py
│   └── superconductivity.py
├── pipelines/          # 多步骤流程
├── io/
│   ├── readers/        # 输出文件解析
│   └── writers/        # 输入文件生成
├── analysis/
│   ├── plotting/       # 绘图
│   ├── superconductivity/  # Tc计算
│   └── thermodynamics/     # 热力学分析
└── utils/
    ├── checkers/       # 任务检查
    └── progress.py     # 进度监控
```

---

## 版本信息

- **当前版本**: 2.0.0
- **CLI版本**: 1.0.0
- **更新日期**: 2025-11-20
- **QE兼容版本**: 6.8+, 7.0+

---

## 相关资源

- [QE官方文档](https://www.quantum-espresso.org/documentation/)
- [QE论坛](https://www.quantum-espresso.org/forum/)
- [PSLibrary赝势](http://pseudopotentials.quantum-espresso.org/)
- [SSSP赝势](https://www.materialscloud.org/discover/sssp)
- [Materials Cloud](https://www.materialscloud.org/)

---

## 技术支持

如遇问题，请：
1. 查看本文档的[常见问题](#常见问题)章节
2. 检查 `--log-level DEBUG` 输出
3. 查阅QE官方文档
4. 在项目Issues中提问

---

**祝你计算顺利！** 🚀
