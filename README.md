# <div align="center"> **my_script**  </div>


# 介绍
该软件是一个计算qe，vasp的小程序。也可以用来产生结构

# 软件架构
qe, vasp, structuregenerator是三个独立的项目，互相不耦合。可以独立开发，使用。


# 安装教程

方法1：
```shell
    python setup.py develop
```
方法2：
```shell
    pip install -e .
```


# <div align="center"> <span style="color:red"> QE篇 </span> </div>
## 基本输入选项。一定要设置的部分：
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

## 具体其它详细的任务模式说明：
**WARNING1 queue存在则会运行，queue不存在则只会产生输入文件和提交任务的脚本文件。**
**WARNING2 如果你使用-j bash, 那么一定注意core设置的不要太大，小心把主节点搞崩溃了。**

###  <span style="color:yellow"> 结构弛豫：
```shell
relax -m mode=relax-vc kpoints_dense="20 20 20" core=48 npool=4  queue=lhy
```
###  <span style="color:yellow"> 自洽
```shell
scf -m mode=scffit kpoints_dense='24 24 24' core=48 npool=4  queue=lhy
```

```shell
scf -m mode=scf kpoints_sparse='12 12 12' core=48 npool=4  queue=lhy
```

###  <span style="color:yellow"> 非自洽计算
```shell
qe_main.py -i relax.out -j bash scf -m mode=nscf kpoints_dense='32 32 32' core=2 queue=local
```

###  <span style="color:yellow"> 不分q点计算声子
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=False queue=lhy core=48 npool=4  queue=lhy el_ph_nsigma=
```

###  <span style="color:yellow"> 分q点计算声子 : split_dyn0模式
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=True core=1 npool=1 queue=local
```
```shell
phono -m mode=split_dyn0 qpoints='6 6 6' core=48 npool=4 queue=local
```

###  <span style="color:yellow"> 分q点计算声子 : split_assignQ模式
```shell
phono -m mode=split_assignQ qpoints='6 6 6' core=1 npool=1 queue=local
```

###  <span style="color:yellow"> 合并声子文件
```shell
phono -m mode=merge core=1 queue=local
```

###  <span style="color:yellow"> 计算力常数
```shell
phono -m mode=q2r qpoints='6 6 6' core=1 npool=1 queue=local
```

###  <span style="color:yellow"> 计算动力学矩阵元
```shell
phono -m mode=matdyn qpoints='6 6 6' core=1 npool=1 queue=local qinserted=50
```

###  <span style="color:yellow"> 计算phonodos, 计算态密度时要用更密的q点网格，这需设置nk1, nk2, nk3   
```shell
dos -m mode=phonodos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 
```

###  <span style="color:yellow"> 计算eletdos(这里计算电子的dos也用qpoints其实非常不合理)
```shell
dos -m mode=eletdos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 
```

###  <span style="color:yellow"> 计算elepdos(这里计算电子的dos也用qpoints其实非常不合理)
```shell
dos -m mode=elepdos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 
```

###  <span style="color:green">使用McAD方法计算超导
####   <span style="color:yellow">不指定最高频率, 将会自动读取最高频率文件
```shell
sc -m mode=McAD core=1 npool=1 queue=local  deguass=0.5 screen_constant=0.1 smearing_method=1 qpoints='6 6 6'
```
####  <span style="color:yellow"> 使用McAD方法超导转变温度指定最高频率
```shell
sc -m mode=McAD core=1 npool=1 queue=local top_freq=80 deguass=0.5 screen_constant=0.1 smearing_method=1 qpoints='6 6 6'
```
这里的参数对应于lamda.in文件中的参数分别是：
```shell
top_freq(最高声子频率)  deguass(展宽宽度取0.12)  smearing_method(展宽方法一般等于1）    
10                               
 0.000000000000000E+00 0.000000000000000E+00 0.000000000000000E+00  1    
 0.000000000000000E+00 0.000000000000000E+00 0.249999999999993E+00  6    
 0.000000000000000E+00 0.000000000000000E+00 -0.499999999999986E+00  3    
 0.000000000000000E+00 0.249999999999993E+00 0.249999999999993E+00  12    
 0.000000000000000E+00 0.249999999999993E+00 -0.499999999999986E+00  12    
 0.000000000000000E+00 -0.499999999999986E+00 -0.499999999999986E+00  3    
 0.249999999999993E+00 0.249999999999993E+00 0.249999999999993E+00  8    
 0.249999999999993E+00 0.249999999999993E+00 -0.499999999999986E+00  12    
 0.249999999999993E+00 -0.499999999999986E+00 -0.499999999999986E+00  6    
 -0.499999999999986E+00 -0.499999999999986E+00 -0.499999999999986E+00  1    
 elph_dir/elph.inp_lambda.1    
 elph_dir/elph.inp_lambda.10    
 elph_dir/elph.inp_lambda.2    
 elph_dir/elph.inp_lambda.3    
 elph_dir/elph.inp_lambda.4    
 elph_dir/elph.inp_lambda.5    
 elph_dir/elph.inp_lambda.6    
 elph_dir/elph.inp_lambda.7    
 elph_dir/elph.inp_lambda.8    
 elph_dir/elph.inp_lambda.9    
screen_constant(库伦屏蔽常数0.1~0.13)
```


### <span style="color:green"> 使用eliashberg方法超导转变温度
### <span style="color:yellow">  使用eliashberg方法超导转变温度, 指定a2F.dos*文件
```shell
sc -m mode=eliashberg core=1 npool=1 queue=local temperature_steps=100 a2F_dos=a2F.dos3 qpoints='6 6 6' screen_constant=0.1
```
####  <span style="color:yellow"> 使用eliashberg方法超导转变温度, 指定使用alpha2F.dat文件中使用哪一列的degauss对应的alpha2F数值。使用degauss_column来指定
(这个方法生成ALPHA2F.OUT可能有问题导致 ELIASHBERG_GAP_T.OUT 中出现NAN。所以更推荐上面那种处理方法。)
```shell
sc -m mode=eliashberg core=1 npool=1 queue=local temperature_steps=100 degauss_column=7 qpoints='6 6 6' screen_constant=0.1
```
NOTE: INPUT文件中只需设置两个参数，
1. 前者是screen_constant，一般取0.10~0.13；
2. 后者是temperature_steps，表示对ntemp个温度点自洽求解Eliashberg方程，得到带隙Δ关于温度T的曲线(该程序首先处理McMillan方程，得到超导临界温度tc作为参考值，然后在温度区间[tc/6, tc*3]中线性插入ntemp个温度点。ntemp一般取40~100即可，也可以更大，建议根据体系差异灵活调控)。


####  <span style="color:yellow"> 获得eliashberg计算得到的超导转变温度
```shell
sc -m mode=eliashberg Tc=output core=1
```

###  <span style="color:yellow">  批量计算
```shell
prepare -m mode=prepare electron_maxstep=1000 core=4 npool=1 queue=local
```
```shell
prepare -m mode=prepare electron_maxstep=1000 core=4 npool=1 queue=local
```


###  <span style="color:yellow"> 如何增加新的功能模块(以增加ele-dos计算的功能模块为例子说明，修改这需要添加哪些内容)

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



# <div align="center"> <span style="color:red"> VASP篇 </span> </div>

## 基本输入选项。一定要设置的部分：
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

## 具体其它详细的任务模式说明：

### 写在最前面，如果你只是希望准备输入文件，即你只是希望准备INCAR,POSCAR,POTCAR,KPOINTS，那么你只需要指定mode与core即可。

### 我会在下面的每一个计算模式下都写好需要该计算模式的最少参数



###  <span style="color:yellow"> 多次结构弛豫到一个离子步收敛 </span>
#### 最简参数
```shell
relax -m mode=rvf core=核数 
```
#### 最繁参数
```shell
relax -m mode=rvf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 三次结构弛豫  </span>
#### 最简参数
```shell
relax -m mode=rv3 core=核数
```
#### 最繁参数
```shell
relax -m mode=rv3 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 单次结构弛豫 </span>
#### 最简参数
```shell
relax -m mode=rv1 core=核数
```
#### 最繁参数
symprec=1e-2 用来解决倒格子对称性和晶格对称性不匹配的问题
```shell
relax -m mode=rv1 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 encut=800 symprec=1e-2
```

###  <span style="color:yellow"> 清理数据, 保留:'POSCAR', 'PPOSCAR', 'POTCAR', 'OUTCAR', 'INCAR*', '*.sh', '*.vasp', '*.slurm'  </span>
```shell
vasp_main.py -w ./ clear -m mode=all
```

###  <span style="color:yellow"> 批量结构弛豫 同时 每个结构多次结构弛豫到一个离子步收敛  </span>
#### 最简参数
```shell
batchrelax -m mode=rvf core=核数
```
#### 最繁参数
```shell
batchrelax -m mode=rvf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 批量结构弛豫 同时 每个结构做三次结构弛豫  </span>
#### 最简参数
```shell
batchrelax -m mode=rv3 core=核数
```
#### 最繁参数
```shell
batchrelax -m mode=rv3 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

### <span style="color:yellow"> 批量结构弛豫 同时 每个结构做三次结构弛豫  </span>
#### 最简参数
```shell
batchrelax -m mode=rv1 core=核数
```
#### 最繁参数
```shell
batchrelax -m mode=rv1 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 有限位移法计算声子谱  </span>
#### 最简参数
```shell
phono -m mode=disp supercell='x x x' core=核数
```
#### 最繁参数
```shell
phono -m supercell='2 2 2' kdensity='36 36 36' mode=disp core=48 ismear=1 encut=800 ediff=1E-08 ediffg=-0.001 queue=lhy
```

###  <span style="color:yellow"> 密度泛函微扰DFPT法计算声子谱  </span>
#### 最简参数
```shell
phono -m mode=dfpt supercell='x x x' core=核数
```
#### 最繁参数
```shell
phono -m supercell='2 2 2' kdensity='36 36 36' mode=dfpt core=48 ismear=1 encut=800 ediff=1e-08 ediffg=-0.001 queue=lhy
```

###  <span style="color:yellow"> 有限位移法计算声子谱——数据处理band  </span>
#### 最简参数
```shell
data -m mode=dispprog supercell='2 2 2' spectrum=True 
```
#### 最繁参数
```shell
data -m mode=dispprog supercell='2 2 2' spectrum=True 
```

###  <span style="color:yellow"> 密度泛函微扰DFPT法计算声子谱——数据处理band  </span>
#### 最简参数
```shell
data -m mode=dfptprog supercell='2 2 2' spectrum=True
```
#### 最繁参数
```shell
data -m mode=dfptprog supercell='2 2 2' spectrum=True
```

###  <span style="color:yellow"> 自洽计算  </span>
#### 最简参数
```shell
properties -m mode=scf core=核数
```
#### 最繁参数
```shell
properties -m mode=scf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800 queue=lhy
```

###  <span style="color:yellow"> 电子态密度计算  </span>
如果指定 "-w 工作路径work_path"， 那么eledos计算时就会在工作路径work_path下寻找一个叫scf的子路径sub_workpath，并将其中的CHGCAR拷贝入eband的子路径sub_workpath
如果没有指定 "-w 工作路径work_path"， 那么默认当前路径是子路径sub_workpath，其母路径就是work_path, 然后重复上述过程， 即：eledos计算时就会在工作路径work_path下寻找一个叫scf的子路径sub_workpath，并将其中的CHGCAR拷贝入eband的子路径sub_workpath

电子态密度计算时，其kspacing需要是scf计算的2倍
#### 最简参数
```shell
properties -m mode=eledos core=核数
```
#### 最繁参数
```shell
properties -m mode=eledos core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.09 encut=800 queue=lhy
```

###  <span style="color:yellow"> 电子能带结构计算  </span>
如果指定 "-w 工作路径work_path"， 那么能带计算时就会在工作路径work_path下寻找一个叫scf的子路径sub_workpath，并将其中的CHGCAR拷贝入eband的子路径sub_workpath
如果没有指定 "-w 工作路径work_path"， 那么默认当前路径是子路径sub_workpath，其母路径就是work_path, 然后重复上述过程， 即：能带计算时就会在工作路径work_path下寻找一个叫scf的子路径sub_workpath，并将其中的CHGCAR拷贝入eband的子路径sub_workpath
#### 最简参数
```shell
properties -m mode=eband core=核数
```
#### 最繁参数
```shell
properties -m mode=eband core=28 ediff=1e-8 ediffg=-0.001 ismear=1 encut=800 queue=lhy
```


# <div align="center"> <span style="color:red"> mytoolkit篇 </span> </div>

## 格式转化
```shell
tool_main.py -i 输入文件名称 -w ./ convert -m dst_format=输出文件名称
```

### POSCAR -> cif, dst_format现在支持的参数为: cif, vasp, struct(wien2k格式)
```shell
tool_main.py -i CaH6.vasp -w ./ convert -m dst_format=CaH6.cif
```



# <div align="center"> <span style="color:red"> 绘制convex hull篇 </span> </div>

## 处理结构优化好的数据

```shell
VaspProcess.py -d ./ -opt
```

## 绘制convex hull

-ebh 10 只关注energy above hull 10meV以内的结构
-hand Lu N H 指定端点
```shell
plot_ternary_convexhull.py -i enthalpy_sorted.csv -hand Lu N H -ebh 10  -save 
```

# <div align="center"> <span style="color:red"> structuregenerator篇 </span> </div>

## 使用splitwps方法产生结构

### 输入文件, 输入文件的名称为splitwps.ini
```python
[splitwps]
spacegroup_number = [[2, 230]]
nameofatoms = ["La", "H"]
numberofatoms = [1, 10]
occupied_number_ranges=[[1,4],[1,10]]
mults_ranges=[[1,2],
              [1,6],]
popsize=100
maxlimit=100
distancematrix=[[2.074, 1.444],
                [1.444, 0.815]]
clathrate_ratio=1.0
hydrogen_content=0.6
max_workers=8
[pso]
numberOflbest = 4
simthreshold = 0.06
fingerprint = "bcm"
lbest = 4
critic = "enthalpy"
maxstep= 50
pso_ltype=["cubic"]
pso_ratio=0.5
```

### 说明：
```shell
max_workers=8 如果指定了，那么就是相当于用8个核进行进程并发
```

### 运行命令:

```shell
generator_main.py -i split.ini -w . method -m mode=splitwps
-i split.ini   输入文件路径
-w .           输出的POSCAR_*的存储路径
mode=splitwps  采用的产生结构的方法
```


## 使用specifywps方法产生结构

### 输入文件, 输入文件的名称一般为specifywps.ini
```python
[splitwps]
spacegroup_number = [[2, 230]]
nameofatoms = ["La", "H"]
numberofatoms = [1, 10]
occupied_number_ranges=[[1,4],[1,10]]
mults_ranges=[[1,2],
              [1,6],]
popsize=100
maxlimit=100
distancematrix=[[2.074, 1.444],
                [1.444, 0.815]]
clathrate_ratio=1.0
hydrogen_content=0.6
max_workers=8
[pso]
numberOflbest = 4
simthreshold = 0.06
fingerprint = "bcm"
lbest = 4
critic = "enthalpy"
maxstep= 50
pso_ltype=["cubic"]
pso_ratio=0.5
```