# my_script

## 介绍
该软件是一个计算qe，vasp的小程序。也可以用来产生结构

## 软件架构
qe, vasp, structuregenerator是三个独立的项目，互相不耦合。可以独立开发，使用。


## 安装教程

方法1：
```shell
    python setup.py develop
```
方法2：
```shell
    pip install -e .
```

## qe使用说明
### 基本输入选项。一定要设置的部分：
```shell
qe_main.py -i 输入文件路径 -w 工作目录 -p 压强 -j 运行方式
```
说明：输入文件路径 和工作路径说明 ：
1. 如果输入文件是relax.out, 那么工作目录会被强行限制为输入文件是relax.out的上一级路径。
2. 如果输入文件是其它路径:
    1. 如果指定了-w ，那么就是在-w指定的路径下创建一个压强值命名的目录，在该压强值命名的目录下开展计算
    2. 如果没有指定-w, 那么就默认所有的计算都在当前指定qe_main.py命令的目录下运行. 并不会额外创建一个压强值命名的目录作为最终工作目录

说明：压强
1. 指定结构优化的压强，单位是GPa.

说明：运行方式
1. bash 代表本地使用bash shell运行。
2. slurm 代表使用slurm脚本运行。
3. pbs 代表使用pbs脚本运行。

说明: 赝势文件最终存放在压强命名的目录的下面

### 具体其它详细的任务模式说明：
**WARNING1 queue存在则会运行，queue不存在则只会产生输入文件和提交任务的脚本文件。**
**WARNING2 如果你使用-j bash, 那么一定注意core设置的不要太大，小心把主节点搞崩溃了。**

结构弛豫：
```shell
relax -m mode=relax-vc kpoints_dense="20 20 20" core=48 npool=4  queue=lhy
```
自洽
```shell
scf -m mode=scffit kpoints_dense='24 24 24' core=48 npool=4  queue=lhy
```
```shell
scf -m mode=scf kpoints_sparse='12 12 12' core=48 npool=4  queue=lhy
```
非自洽计算
```shell
qe_main.py -i relax.out -j bash scf -m mode=nscf kpoints_dense='32 32 32' core=2 queue=local
```

不分q点计算声子
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=False queue=lhy core=48 npool=4  queue=lhy el_ph_nsigma=
```
分q点计算声子 : split_dyn0模式
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=True core=1 npool=1 queue=local
```
```shell
phono -m mode=split_dyn0 qpoints='6 6 6' core=48 npool=4 queue=local
```
分q点计算声子 : split_assignQ模式
```shell
phono -m mode=split_assignQ qpoints='6 6 6' core=1 npool=1 queue=local
```

合并声子文件
```shell
phono -m mode=merge core=1 queue=local
```

计算力常数
```shell
phono -m mode=q2r qpoints='6 6 6' core=1 npool=1 queue=local
```
计算动力学矩阵元
```shell
phono -m mode=matdyn qpoints='6 6 6' core=1 npool=1 queue=local qinserted=50
```
计算态密度, 计算态密度时要用更密的q点网格，这需设置nk1, nk2, nk3   
```shell
dos -m mode=dos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 
```

计算超导
不指定最高频率
```shell
sc -m mode=McAD core=1 npool=1 queue=local deguass=0.5 screen_constant=0.1 smearing_method=1 qpoints='6 6 6'
```
指定最高频率
```shell
sc -m mode=McAD core=1 npool=1 queue=local top_freq=80 deguass=0.5 screen_constant=0.1 smearing_method=1 qpoints='6 6 6'
```
```shell
sc -m mode=eliashberg core=1 npool=1 queue=local temperature_points=10000 a2F_dos=a2F.dos3 qpoints='6 6 6'
```

批量计算
```shell
prepare -m mode="relax-vc scffit scf nosplit" dyn0_flag=True qpoints='6 6 6' core=4 npool=1 queue=local
```
```shell
prepare -m mode="relax-vc scffit scf        " dyn0_flag=True qpoints='6 6 6' core=4 npool=1 queue=local
```
### 如何增加新的功能模块(以增加ele-dos计算的功能模块为例子说明，修改这需要添加哪些内容)

#### 第1步：在qe_writeinput.py中, class qe_writeinput中增加写电子态密度计算.in文件的实例方法
#### 第2步：在qe_writeinput.py中, class qe_writeinput的dos的init_from_dosinput的类方法中，新增需要用到的新的输入参数
```python
DeltaE=other_class.DeltaE,
emin=other_class.emin,
emax=other_class.emax,
```
#### 第3步: 在qedos_inputpara.py中补充关于新的输入参数的默认变量设置class qedos_inputpara(qe_inputpara)
```python
# 电子态密度设置
if not hasattr(self, "DeltaE"):
    self.DeltaE = 0.01
    logger.warning("You didn't set `DeltaE`, the program will use default value: DeltaE=0.01")

# 电子态密度设置
if not hasattr(self, "emin"):
    self.emin = -10
    logger.warning("You didn't set `emin`, the program will use default value: emin=-10")

# 电子态密度设置
if not hasattr(self, "emax"):
    self.emax = 30
    logger.warning("You didn't set `emax`, the program will use default value: emax=30 ")
```
#### 第4步: 

## vasp使用说明

### 基本输入选项。一定要设置的部分：
```shell
vasp_main.py -i 输入文件路径 -w 工作目录 -p 压强 -j 运行方式
```
说明：输入文件路径：
1. 指定输入的POSCAR

说明：工作目录
1. 如果指定工作目录，那么工作目录是-w指定的目录+-p指定的压强组成的路径。比如: -w test -p 200。那么最终路径就是/test/200.0。所有的文件都会在这个路径下运行。
2. 如果没有指定工作路径，那么工作路径就是当前路径，不会产生压强目录。

说明：压强
1. 指定结构优化的压强，单位是GPa.

说明：运行方式
1. bash 代表本地使用bash shell运行。
2. slurm 代表使用slurm脚本运行。
3. pbs 代表使用pbs脚本运行。

### 具体其它详细的任务模式说明：

清理数据, 保留:'POSCAR', 'PPOSCAR', 'POTCAR', 'OUTCAR', 'INCAR*', '*.sh', '*.vasp', '*.slurm'
```shell
vasp_main.py -w ./ clear -m mode=all
```

批量结构弛豫
```shell
vasp_main.py -i ./ -w ./ -p 200 -j slurm batchrelax -m mode=rv3 core=1 
```

```shell
vasp_main.py -i ./POSCAR -j bash phono -m supercell='2 2 2' kpoints='2 2 2' mode=disp core=1 queue=local
```

disp声子计算
```shell
vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm phono -m supercell='2 2 2' kpoints='40 40 40' mode=disp  
```

dfpt声子计算
```shell
vasp_main.py -i ./test/POSCAR -w ./test/ -j slurm phono -m supercell='2 2 2' kpoints='40 40 40' mode=dfpt
```

## structuregenerator使用说明


