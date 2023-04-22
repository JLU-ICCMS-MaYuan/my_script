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
1. 如果输入文件是relax.out或者scffit.out或者scf.out, 那么工作目录会被强行限制为输入文件所在目录
2. 如果输入文件是类型例如，POSCAR或者其它scf.in或者cif文件:
    1. 如果指定了-w ，那么就是在-w指定的路径下开展计算
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

结构优化时，对于体系的对称性非常敏感，也就是QE自身不会寻找体系的对称性，只能依靠手动输入体系的对称性，这点与VASP相比是非常欠缺的，因为VASP是可以自己寻找对称性的
```shell
relax -m mode=relax-vc kpoints_dense="20 20 20" conv_thr=1.0d-8 core=48 npool=4  queue=lhy
```
###  <span style="color:yellow"> 自洽
#### 重要参数的说明
1. mixing_beta = 0.3 一般用0.3，算的比较快
2. degauss = 0.05, 一般用0.02~0.05，值越大，算出来的声子谱越容易稳定。特别是当计算出来的声子谱有一个非GAMMA点虚频一点点时，可以通过调整degauss来消除虚频

**所以我采用的默认参数mixing_beta=0.3, degauss=0.05, 如果你自己需要高精度，自己调整**
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
### <span style="color:yellow"> 声子计算
#### 重要参数的说明
1. alpha_mix(1)=0.3 一般用0.3， 算的比较快，值越大算的越慢。

**所以我采用的默认参数alpha_mix(1)=0.3, 如果你自己需要高精度，自己调整**

####  不分q点计算声子
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=False queue=lhy core=48 npool=4  queue=lhy el_ph_nsigma=
```

####  <span style="color:green"> 分q点计算声子 : split_dyn0模式
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=True core=1 npool=1 queue=local
```
```shell
phono -m mode=split_dyn0 qpoints='6 6 6' core=48 npool=4 queue=local
```

####  <span style="color:green"> 分q点计算声子 : split_assignQ模式
```shell
phono -m mode=split_assignQ qpoints='6 6 6' core=1 npool=1 queue=local
```

####  <span style="color:green"> 合并声子文件
```shell
phono -m mode=merge core=1 queue=local
```

####  <span style="color:green"> 计算力常数
```shell
phono -m mode=q2r qpoints='6 6 6' core=1 npool=1 queue=local
```

####  <span style="color:green"> 计算动力学矩阵元，获得高对称路径下的声子振动频率
```shell
phono -m mode=matdyn qpoints='6 6 6' core=1 npool=1 queue=local qinserted=50
```

####  <span style="color:green"> 处理数据获得声子谱
**注意！！！ 计算完高对称路径下的声子振动频率后，一定要记得处理数据获得声子谱。输出一个文件 qp_freq_width.csv 可以用来绘制声子谱以及每个振动频率对应的线宽**
这里是处理声子谱的命令
```shell
phono -m mode=phonobanddata core=1
```
**这里是获得高对称点投影路径的命令， 可以用来在origin中绘制高对称点的垂直参考线**
```shell
phono -m mode=hspp core=1
```

**不然如果后续直接计算声子态密度phonodos的话，`体系名称`.freq文件, `体系名称`.freq.gq文件中声子的q点数目就和gam.lines文件中q点数目对不上了。**
这里展示计算完`动力学矩阵元，获得高对称路径下的声子振动频率`后的 `体系名称`.freq文件(因为`体系名称`.freq文件和`体系名称`.freq.gq文件类似，这里就不再展示）, gam.lines文件

我先给出我的matdyn.in文件: 其中有6个高对称的q点，每2个q点之间插入20个q点，所以总的q点个数是：20*5+1=101， 101个q点。
```shell
&input    
 asr = 'simple',    
 amass(1)=138.90547 ,    
 amass(2)=88.90585 ,    
 amass(3)=1.00794 ,    
  flfrc = 'La1Y3H40.fc',    
  flfrq='La1Y3H40.freq',    
  la2F=.true.,    
  dos=.false.,    
  q_in_band_form=.true.,    
  q_in_cryst_coord=.true.,    
/    
6    
 0               0               0               20    
 0.0             0.5             0.0             20    
 0.5             0.5             0.0             20    
 0               0               0               20    
 0.5             0.5             0.5             20    
 0.0             0.5             0.0             20  
```
`体系名称`.freq文件
```shell
&plot nbnd= 132, nks= 101 / # nbnd代表有132个声子振动模式(反应在声子谱上有132条带)，nks代表计算了101个q点
           0.000000  0.000000  0.000000 # 第1个q点的坐标， 下面的这些数字正好是132个，表示在这个q点下计算了132个振动模式的振动频率分别是多少。
-0.0663   -0.0663   -0.0663  188.7484  188.7484  188.7484
...       ...       ...      ...       ...       ...
...       ...       ...      ...       ...       ...
2044.5579 2044.5579 2226.2463 2226.2463 2253.1870 2327.3298
           0.000000  0.025000  0.000000 # 第2个q点的坐标
6.9385    6.9385   25.9481  188.0445  188.0445  189.0030
...       ...       ...      ...       ...       ...
...       ...       ...      ...       ...       ...
2044.9873 2045.0136 2225.8874 2226.5873 2253.9177 2328.7849
...
...
            0.000000  0.500000  0.000000 # 第101个q点的坐标
48.9448   48.9448  140.2440  140.2440  227.4149  227.4149
...       ...       ...      ...       ...       ...
...       ...       ...      ...       ...       ...
1996.1023 2015.9340 2157.0716 2191.5706 2244.4918 2301.3597 
```
gam.lines文件
```shell
   
  Gamma lines for all modes [THz] 
   
 Broadening   0.0050 # degauss展宽，第1个值是起始展宽值，
                     # 取决于你的ph.in中el_ph_sigma=0.005和el_ph_nsigma=10的设置
       1  # 代表第1个q点, 
  # 下面的数字表示在该q点下，每个振动频率的展宽值，因为每个q点有132个振动模式，每个振动模式有一个展宽，所有总共有132个值
  0.0000  0.0000  0.0000  0.0075  0.0075  0.0075  0.0152  0.0152  0.0152 
  ...     ...     ...     ...     ...     ...     ...     ...     ...     
  ...     ...     ...     ...     ...     ...     ...     ...     ... 
  0.9044  0.9044  1.7802  1.7802  6.8412 11.6098
       101
  0.0479  0.0479  0.0049  0.0049  0.0027  0.0027  0.0259  0.0046  0.0046
  ...     ...     ...     ...     ...     ...     ...     ...     ...     
  ...     ...     ...     ...     ...     ...     ...     ...     ... 
  0.2878  1.5340  0.7457  2.5516  3.0478  2.3681
 Broadening   0.0100 # degauss展宽，第2个展宽
       1
       ...
       101
       ...
 Broadening   0.0500
       1
       ...
       101
       ...
```

### 声子态密度计算
####  <span style="color:green"> 计算phonodos, 计算态密度时要用更密的q点网格，这需设置nk1, nk2, nk3   
```shell
phono -m mode=phonodos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 
```

####  <span style="color:green"> 处理出投影到每种元素的phonodos，并且会输出一个文件phdos_proj2eles.csv用来绘制投影态密度
```shell
phono -m mode=phonodosdata core=1 queue=local
```

###  <span style="color:yellow"> 电子态密度计算
####  <span style="color:green">**总DOS: eletdos**</span>
```shell
eletron -m mode=eletdos core=1 npool=1 queue=local kpoints_dense='8 8 8' 
```

####  <span style="color:green">**投影DOS: elepdos**</span>
```shell
eletron -m mode=elepdos core=1 npool=1 queue=local kpoints_dense='8 8 8' 
```

###  <span style="color:yellow"> 电子能带结构计算
####  <span style="color:green">**获得eleband数据**</span>
```shell
eletron -m mode=eleband core=1 npool=1 queue=local kinserted=200 nbnd=500
```

####  <span style="color:green">**处理eleband数据获得可以origin绘图的数据**</span>(这里计算电子的dos也用qpoints其实非常不合理)
```shell
eletron -m mode=elebanddata core=1 queue=local
```

###  <span style="color:yellow"> 同时 电子能带结构计算 电子态密度计算
```shell
eletron -m mode=eleproperties core=48 npool=4 queue=local kinserted=200 nbnd=500  kpoints_dense='8 8 8' 
```
####  <span style="color:green">**处理eleband数据获得可以origin绘图的数据**</span>(这里计算电子的dos也用qpoints其实非常不合理)
```shell
eletron -m mode=elebanddata core=1 queue=local
```


###  <span style="color:yellow"> 计算超导
处理电荷屏蔽常数为0.1和0.13,得到 $\lambda$ 和 $\omega_{log}$ 并且输出一个文件：w_alpha2f_lambda.csv 可以用来绘制alpha-lambda的函数图像

```shell
sc -m mode=Tc core=1 npool=1 queue=local temperature_steps=100 qpoints='4 4 4' a2fdos=True alpha2fdat=False broaden=0.5 smearing_method=1 gaussid=3 gauss=0.015  top_freq=80
```

```shell
sc -m mode=Tc core=1 npool=1 queue=local temperature_steps=100 qpoints='4 4 4' a2fdos=True alpha2fdat=False broaden=0.5 smearing_method=1
```
指定最高频率
```shell
top_freq=80
```
**当然了，a2fdos， alpha2fdat 两个参数不设置也可以，默认使用的是从alpha2fdat中读取数据，因为很多时候，你不会计算phonodos，所以你也没有a2F.dos*这些文件。**
使用eliashberg方法超导转变温度, 指定读取a2F.dos*文件
```shell
a2fdos=True alpha2fdat=False
```
使用eliashberg方法超导转变温度, 指定读取alpha2F.dat文件中使用哪一列的degauss对应的alpha2F数值。使用alpha2fdat来指定
**(这个方法生成ALPHA2F.OUT可能有问题导致 ELIASHBERG_GAP_T.OUT 中出现NAN。所以更推荐a2fdos=True那种处理方法。)**
```shell
a2fdos=False alpha2fdat=True
```

用到的公式：

$$\lambda = 2 \int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega$$

$$\omega_{log} = exp[\frac{2}{\lambda}\int_{0}^{\infty} \frac{d\omega}{\omega} \alpha^2 F(\omega) ln(\omega)]$$

####  <span style="color:green"> **千万注意：alpha2F.dat中频率的单位是THz, 但是freq.gp文件中的频率的单位是cm-1，如果想在一张图中把两个数据放在一起对比需要一个单位转化**
$$c=299792458  m/s$$
$$\lambda^{-1} = \frac{\nu}{c} = \frac{1Thz}{299792458  m/s} =  \frac{10^{12}Hz}{299792458\times 10^{2} cm \cdot Hz} = 33.3564095198152 cm^{-1} $$
即：
$$1Thz \Leftrightarrow 33.3564 cm^{-1}$$

w_alpha2f_lambda.csv 中的第一列是频率，其单位是经过转化的$cm^{-1}$


####  <span style="color:green"> McAD方法计算Tc需要lamda.in文件中的参数分别是：
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

####  <span style="color:green"> Eliashberg方法计算Tc需要INPUT文件中只需设置两个参数，
1. 前者是screen_constant，一般取0.10~0.13；
2. 后者是temperature_steps，表示对ntemp个温度点自洽求解Eliashberg方程，得到带隙Δ关于温度T的曲线(该程序首先处理McMillan方程，得到超导临界温度tc作为参考值，然后在温度区间[tc/6, tc*3]中线性插入ntemp个温度点。ntemp一般取40~100即可，也可以更大，建议根据体系差异灵活调控)。



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
1. 如果指定工作目录，那么工作目录是-w指定的目录
2. 如果没有指定-w, 那么就默认所有的计算都在当前指定vasp_main.py命令的目录下运行. 并不会额外创建一个压强值命名的目录作为最终工作目录

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
-w 你的路径/scf properties -m mode=scf core=核数
```
#### 最繁参数
```shell
-w 你的路径/scf properties -m mode=scf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800 queue=lhy
```

###  <span style="color:yellow"> 电子态密度计算  </span>

如果指定 "-w 工作路径work_path"， 那么eledos计算时就会在工作路径work_path的母路径下寻找一个叫scf的目录，并将scf其中的CHGCAR拷贝入指定的路径work_path下; 所以在指定work_path时，尽量写出来目录eledos

如果没有指定 "-w 工作路径work_path"， 那么eledos计算时就会在当前路径的母路径下寻找一个叫scf的目录，并将scf其中的CHGCAR拷贝入当前路径work_path下。

**电子态密度计算时，INCAR中有几个参数需要格外注意，这里我将这些参数给出，请你使用该脚本产生eledos计算的INCAR后稍作检查**
```shell
ISTART = 1    # 读取WAVECAR，如果WAVECAR不存在，就自己构造
ICHARG = 11   # 读取CHGCAR， 进行非自洽计算
#EMIN = -10   # 该脚本中默认将其注释了，EMIN指定了DOS评估的能量范围的下限。
#EMAX =  10   # 该脚本中默认将其注释了，EMAX指定了DOS评估的能量范围的上限。
```

电子态密度计算时，其kspacing需要是scf计算的2倍
#### 最简参数
```shell
-w 你的路径/eledos properties -m mode=eledos core=核数
```
#### 最繁参数
```shell
-w 你的路径/eledos properties -m mode=eledos core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.09 encut=800 queue=lhy
```

###  <span style="color:yellow"> 电子能带结构计算  </span>

如果指定 "-w 工作路径work_path"， 那么eband计算时就会在工作路径work_path的母路径下寻找一个叫scf的目录，并将scf其中的CHGCAR拷贝入指定的路径work_path下; 所以在指定work_path时，尽量写出来目录eband。

如果没有指定 "-w 工作路径work_path"， 那么eband计算时就会在当前路径的母路径下寻找一个叫scf的目录，并将scf其中的CHGCAR拷贝入当前路径work_path下。

**电子态密度计算时，INCAR中有几个参数需要格外注意，这里我将这些参数给出，请你使用该脚本产生eledos计算的INCAR后稍作检查**
```shell
ISTART = 0    # 读取WAVECAR，如果WAVECAR不存在，就自己构造
ICHARG = 11   # 读取CHGCAR， 进行非自洽计算
```

#### 最简参数
```shell
-w 你的路径/eband properties -m mode=eband core=核数
```
#### 最繁参数
```shell
-w 你的路径/eband properties -m mode=eband core=28 ediff=1e-8 ediffg=-0.001 ismear=1 encut=800 queue=lhy
```
#### 获得投影到一维路径的高对称路径点
```shell
vasp_main.py -i CONTCAR -j bash data -m mode=hspp core=1
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
##安装注意事项：
由于用到了CALYPSO内的部分Fortran代码，所以在忻宇师弟的帮助下，将Fortran代码进行打包，编译为Python可以调用的动态链接库（几个`.so文件`）。然后将这些动态链接库文件复制到指定的目录中。
### 第一步编译

```shell
cd structuregenerator/psolib/fortran_src/LegacyFingerPrint/
make
# 编译好的文件存储在 structuregenerator/psolib/fortran_src/sym3dgenerator/LIB
ls LIB
_f90sym3dgenerator.cpython-39-x86_64-linux-gnu.so f90sym3dgenerator.cpython-39-x86_64-linux-gnu.so
cd structuregenerator/psolib/fortran_src/sym3dgenerator
make
ls LIB
_f90LegacyFingerPrint.cpython-39-x86_64-linux-gnu.so f90LegacyFingerPrint.cpython-39-x86_64-linux-gnu.so
```
### 第二步复制LIB目录中`.so文件`到指定目录
```shell
cd structuregenerator
# 这里可能需要你自己在psolib目录中创建一个新子目录：lib_so
mkdir -p psolib/lib_so
cd psolib/lib_so
cp ../fortran_src/LegacyFingerPrint/LIB/*.so  .
cp ../fortran_src/sym3dgenerator/LIB/*.so  .
```


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
remain_H_ratio_UPPERSTD = 0.44
fraction_of_hydrogen_volume_LOWERSTD = 0.3
shr_num_avg_LOWERSTD = 1.6
cage_regularity_avg_UPPERSTD = 0.05
h2h_network_regularity_avg_UPPERSTD = 0.15
h2h_1nearest_LOWERSTD = 0.9
nonh2nonh_1nearest_LOWERSTD = 2.2
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