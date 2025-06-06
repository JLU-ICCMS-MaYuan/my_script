# <div align="center"> **my_script**  </div>

# 网络热帖:
1. https://zhuanlan.zhihu.com/p/427794442 如何安装qe
2. https://zhuanlan.zhihu.com/p/547184583 如何安装vasp

# 介绍
该软件是一个计算qe，vasp的小程序。也可以用来产生结构

# 软件架构
qe, vasp, structuregenerator是三个独立的项目，互相不耦合。可以独立开发，使用。


# 安装教程

## 第一步：配置my_scriptrc.py 
```python
qebin_path # qe的bin目录
epwbin_path # epw的bin目录
qe_pseudopotential_dir # qe的赝势目录
eliashberg_x_path # eliashberg_x的文件路径
vaspstd_path # vasp_std的文件路径
vaspgam_path # vasp_gam的文件路径
potcar_dir_path # vasp赝势路径

bashtitle中的环境变量, slurmtitle, pbstitle, lsftitle中挑出来你要用的机器支持的系统，然后把对应的环境变量修改好。

# 然后执行
cp my_scriptrc.py  ~/.my_scriptrc.py
python ~/.my_scriptrc.py
```

## 第二步：在my_script目录中执行write_script_path.py，必须在my_script目录中执行！！！
```shell
python write_script_path.py
```

## 第三步：安装软件
方法1（不推荐）：
```shell
python setup.py develop
```
方法2（推荐）：
```shell
# 在有网络的机器上安装该软件
pip install -e .

# 或者在没有网络的机器上安装该软件
# 先在其它有网络的机器上安装好该软件，并将相应的conda环境打包到无网络的机器上，然后执行下面的代码即可。
pip install -e .  --no-build-isolation --no-index --find-links=./
# `--no-index`：不要从PyPI下载包。
#   这个选项告诉pip在安装包时不要使用Python包索引(PyPI)。也就是说，pip不会从PyPI下载包，而是只从指定的本地或网络位置安装包。这通常用于安装那些不在PyPI上发布的包，或者需要从特定位置安装包的情况。
# `--find-links ./`：从指定的本地路径`./`查找包。
#   这个选项后面跟的是一个或多个URL或文件系统路径，pip会从这些指定的位置查找包。这可以是本地文件路径，也可以是网络URL。当你使用`--find-links`时，pip会首先检查这些链接，而不是默认的PyPI索引。

# 你也可以在有网络的机器上执行下面的代码
pip freeze > requirements.txt
# 然后删掉requirements.txt里面关于网址的部分，不删除会报错
# 
pip download --platform anylinux_x86_64 --no-deps on -d pip_packages/ -r requirements.txt
# 一个关于pip简单的教程：https://blog.csdn.net/Leon_Jinhai_Sun/article/details/139171646

```

有时候你安装的软件会有问题，比如：`ModuleNotFoundError: No module named 'qe'`, 这个时候你可以执行下面的代码：
```shell
python -c "import sys; print(sys.path)" 

['', '/data/software/intel/oneapi2021/advisor/2021.1.1/pythonapi', '/data/home/liuhanyu/workplace/luojie/soft/miniconda3/envs/my_scripts/lib/python39.zip', '/data/home/liuhanyu/workplace/luojie/soft/miniconda3/envs/my_scripts/lib/python3.9', '/data/home/liuhanyu/workplace/luojie/soft/miniconda3/envs/my_scripts/lib/python3.9/lib-dynload', '/data/home/liuhanyu/workplace/luojie/soft/miniconda3/envs/my_scripts/lib/python3.9/site-packages']

# 检查这个列表中有没有包含my_scripts的路径，如果没有，那么你需要在你的~/.bashrc文件中添加下面的代码：
export PYTHONPATH="/data/home/liuhanyu/code/my_script:$PYTHONPATH"
```
**但是本质上这个问题的出现是由于`my_script.egg-info/top_level.txt`中没有相应的`qe`为名称的目录名
在`pyproject.toml`中添加`[tool.setuptools] packages = ["qe", "vasp", "epw", "structuregenerator"]`即可解决。此时再打开`my_script.egg-info/top_level.txt`
就会发现，出现了qe, vasp, epw, structuregenerator四个目录名**


# <div align="center"> <span style="color:red"> QE篇 </span> </div>

## 基本输入选项。一定要设置的部分：
```shell
qe_main.py -i 输入文件路径 -w 工作目录 -p 压强 -j 运行方式 -l 输出日志文件等级
```
说明：输入文件路径 和工作路径说明 ：
1. -i, 输入文件任意类型都可以，读入后会拷贝到当前目录以备份
2. -w, 如果指定了-w ，那么就是在-w指定的路径下开展计算. 如果没有指定-w, 那么就默认所有的计算都在当前指定qe_main.py命令的目录下运行. 并不会额外创建一个压强值命名的目录作为最终工作目录
    说明: 赝势文件最终存放在-w命名的目录的下面
3. -p, 只有结构优化该参数起作用，否则其他模式一律无用。单位：GPa
4. -j, 支持的参数slurm, pbs, lsf, bash
    说明：运行方式:
    1. bash 代表本地使用bash shell运行。
    2. slurm 代表使用slurm脚本运行。
    3. pbs 代表使用pbs脚本运行。
5. -l, 支持的参数debug, info, warning, error, critical, 默认为info




## 具体其它详细的任务模式说明：
<span style="color:red"> **WARNING1 queue存在则会运行，queue不存在则只会产生输入文件和提交任务的脚本文件。**

<span style="color:red"> **WARNING2 如果你使用-j bash, 那么一定注意core设置的不要太大，小心把主节点搞崩溃了。**

<span style="color:red"> **WARNING3 在设置q点之前一定要非常谨慎，你设置的q点的密度，k点的密度必须与体系的对称性相匹配，必须匹配，所以建议你用vaspkit生成一下KPOINTS检查一下k点密度，以及用kmesh.py这个小脚本也生成一下k点密度，两者进行对比，最终确定你要用哪个k点密度。这个值最终要放在文章里，所以一定要非常谨慎！！！**

###  <span style="color:yellow"> 结构弛豫：

结构优化时，对于体系的对称性非常敏感，也就是QE自身不会寻找体系的对称性，只能依靠手动输入体系的对称性，这点与VASP相比是非常欠缺的，因为VASP是可以自己寻找对称性的
```shell
relax -m mode=relax-vc kpoints_dense="20 20 20" conv_thr=1.0d-8 execmd='mpirun -np 48' npool=4  queue=lhy
```

###  <span style="color:yellow"> 自洽
#### 重要参数的说明
1. mixing_beta = 0.3 一般用0.3，算的比较快
2. degauss = 0.05, 一般用0.02~0.05，值越大，算出来的声子谱越容易稳定。特别是当计算出来的声子谱有一个非GAMMA点虚频一点点时，可以通过调整degauss来消除虚频

**所以我采用的默认参数mixing_beta=0.3, degauss=0.05, 如果你自己需要高精度，自己调整**
```shell
scf -m mode=scffit kpoints_dense='24 24 24' execmd='mpirun -np 48' npool=4  queue=lhy
```

```shell
scf -m mode=scf kpoints_sparse='12 12 12' execmd='mpirun -np 48' npool=4  queue=lhy
```

###  <span style="color:yellow"> 非自洽计算
```shell
qe_main.py -i relax.out -j bash scf -m mode=nscf kpoints_dense='32 32 32' execmd='mpirun -np 48' queue=local
```

### <span style="color:yellow"> 声子计算

####  <span style="color:green"> 不分q点计算声子
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=False queue=lhy execmd='mpirun -np 48' npool=4  queue=lhy el_ph_nsigma=10 
```
<span style="color:green"> **提示：如果你的工作目录中已经存在了dyn0文件，那么可以不用输入qpoints这个参数，程序会自动读入dyn0中的第一行以获得qpoints**

####  <span style="color:green"> 分q点计算声子 : split_dyn0模式
```shell
phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=True execmd='' npool=1 queue=local
```
```shell
phono -m mode=split_dyn0 qpoints='6 6 6' execmd='mpirun -np 48' npool=4 queue=local
```

####  <span style="color:green"> 分q点计算声子 : split_assignQ模式
```shell
phono -m mode=split_assignQ qpoints='6 6 6' execmd='mpirun -np 48' npool=1 queue=local
```

####  <span style="color:green"> 合并声子文件
```shell
phono -m mode=merge execmd='mpirun -np 48' queue=local
```

####  <span style="color:green"> 合并某一个q点下的所有不可约表示
首先讲解一下某一个q点目录下的tmp目录的文件结构
```shell
# 文件
[prefix].a2Fsave # 非常重要，这个文件也要拷贝到总的tmp目录中。它产生自scffit.in计算中，因为在scffit.in中开启了la2f=.true.这个开关，所以会在tmp目录中产生该文件。在完成tmp/_ph0/[prefix].phsave/中所有动力学矩阵和电声耦合矩阵的读取后，将会读取[prefix].a2Fsave文件进行电声耦合计算。
[prefix].wfc*
[prefix].xml
# 目录
[prefix].save # 该目录下存储着赝势，波函函数文件等
    *.UPF
    wfc*.dat
    charge-density.dat
    data-file-schema.xml
    paw.txt 
_ph0 # 该目录存储着dv*文件， wfc文件， [prefix].phsave目录， [prefix].save目录
    [prefix].[prefix].dv*
    [prefix].[prefix].dv_paw*
    [prefix].wcf*
    [prefix].phsave # 该目录下存储着control_ph.xml，dynmat.1.*.xml，elph.1.*.xml, patterns.1.xml, status_run.xml
        dynmat.1.*.xml
        elph.1.*.xml
        # 一般来说dynmat.1.*.xml的文件数等于elph.1.*.xml的文件数等于不可约表示的数目，
        # 需要通过检查dynmat.1.*.xml中的</Root>关键字来
        # 判断是否左右的不可约表示的动力学矩阵都计算完毕。
        patterns.1.xml
        status_run.xml
        control_ph.xml
        # status_run.xml 和 control_ph.xml 这两个文件一定不可以删掉
        # 这两个文件如果删掉了，即使ph.in中recover=.true.,也会从头计算声子，无法续算。
    [prefix].save 
        # 如果你在中间某个不可约表示计算时出错，中断了计算，
        # 后续采用指定起始不可约表示的方法（start_irr=***)完成了该q点所有表示的计算，
        # 那么这个目录中将只有charge-density.dat  data-file-schema.xml两个文件。
        # 但是后来，我检查了其它体系分q点计算的tmp/_ph0/[prefix].save目录
        # 发现 tmp/_ph0/[prefix].save目录 中也只有charge-density.dat  data-file-schema.xml两个文件
        charge-density.dat  
        data-file-schema.xml
        # 如果你的在合并不可约表示、准备处理出dyn*文件和elph_dir目录时，报了一些错，那么有可能是tmp/_ph0/[prefix].save目录中缺少文件，只需要把tmp/[prefix].save中所有的文件拷贝到tmp/_ph0/[prefix].save中即可。
```

下面讲解一下合并某一个q点的所有不可约表示的步骤：

首先给出合并所有不可约表示所需要的全部文件：###表示该文件是合并所需要的文件

```shell
[prefix].a2Fsave           ###
[prefix].wfc*
[prefix].xml
[prefix].save              ###
    *.UPF                  ###
    wfc*.dat               ###
    charge-density.dat     ###
    data-file-schema.xml   ###
    paw.txt                ###
_ph0                       ###
    [prefix].[prefix].dv*
    [prefix].[prefix].dv_paw*
    [prefix].wcf*
    [prefix].phsave        ###
        dynmat.1.*.xml     ###
        elph.1.*.xml       ###
        patterns.1.xml     ###
        status_run.xml     ###
        control_ph.xml     ###
    [prefix].save          ###
        charge-density.dat ###
        data-file-schema.xml
```

0. 第0步: 修改split_ph.in文件 以及 提交作业的脚本

```shell
删除  start_irr=***, 
添加  recover=.true.,
注释所有关于自洽计算的内容
```

1. 第1步: 检查文件
```shell
cd tmp/_ph0/La2Ce2Y2H54.phsave/
grep "</Root>" dynmat.1.*.xml | wc -l
grep "</Root>" elph.1.*.xml | wc -l 
grep "</root>" elph.1.*.xml | wc -l 

# 最好这两个grep的结果是一样的，但是有时候不一样，有时候会有一个dynmat.1.0.xml文件。

sbatch s5_PhSplitDyn0.sh 
# 尝试第1次合并该q点的所有表示
```

2. 第2步：如果第1步没有问题，那么就不必执行第2步。如果第1步有问题，比如爆出这样的错误
```shell
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine phq_readin (1):
     Electron-phonon only for metals
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     stopping ...
     read_file: Wavefunctions in collected format not available

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     stopping ...
```

那么只需要把tmp/[prefix].save中所有的文件拷贝到tmp/_ph0/[prefix].save中即可
```shell
cd tmp
cp La2Ce2Y2H54.save/* _ph0/La2Ce2Y2H54.save/

sbatch s5_PhSplitDyn0.sh 
# 尝试第2次合并该q点的所有表示
```

3. 第3步：如果还是无法合并所有不可约q点，并且在ph.out中报出这样的警告和错误后停止运行ph.x。这多半是因为你编译的qe有问题。再重新编译qe后虽然这样的错误还会报出，但是不会中止程序运行。
```shell
warning: file closed at level 1 with tag Root open
 end of file reached, tag PARTIAL_EL_PHON not found
 end of file reached, tag NUMBER_OF_K not found
 end of file reached, tag NUMBER_OF_BANDS not found
 end of file reached, tag K_POINT.1 not found
 end of file reached, tag K_POINT.2 not found
 end of file reached, tag K_POINT.3 not found
xmlr_closetag: severe error, closing tag that was never opened
```
4. 第4步：如果还是无法合并所有不可约q点，程序停在了electron-phonon interaction这里，虽然没有报出CRASH，但是在log.err(slurm作业系统提供的错误报告)文件中提示无法读取[prefix].a2Fsave文件。这是因为你缺少[prefix].a2Fsave文件，从某一个不可约表示计算目录的tmp文件中拷贝一个过来就好了。

```shell
     freq ( 176- 176) =       2060.7  [cm-1]   --> A_g
     freq ( 177- 177) =       2096.1  [cm-1]   --> A_u
     freq ( 178- 178) =       2102.2  [cm-1]   --> A_g
     freq ( 179- 179) =       2177.1  [cm-1]   --> A_g
     freq ( 180- 180) =       2183.9  [cm-1]   --> A_u
     electron-phonon interaction  ...
```


###  <span style="color:yellow"> 计算力常数
```shell
phono -m mode=q2r qpoints='6 6 6' execmd='' npool=1 queue=local
```

###  <span style="color:yellow"> 声子态密度计算
####  <span style="color:green"> 计算phonodos, 计算态密度时要用更密的q点网格, 通过设置qpoints获得更密集得q点网格
```shell
phono -m mode=phonodos execmd='mpirun -np 48' npool=1 queue=local qpoints='8 8 8' ndos=500 
```


####  <span style="color:green"> 处理出投影到每种元素的phonodos，并且会输出一个文件phdos_proj2eles.csv用来绘制投影态密度
```shell
phono -m mode=phonodosdata execmd='' queue=local
```

###  <span style="color:yellow"> 考虑振动熵后的吉布斯自由能的计算

```shell
phono -m mode=gibbsvb execmd='mpirun -np 48' queue=local
```
运行完该命令后会得到三个文件分别为： thermodynamics_from_dyn1.csv，thermodynamics_from_dyns.csv，thermodynamics_from_phtdos.csv。他们的文件格式都一样，第一列是温度，第二列是每个结构的自由能，单位是eV，第三列是平均到每原子的自由能，单位是meV。第一行是零点振动能（Zero Point ENERGY）

**thermodynamics_from_dyn1.csv 代表自由能是通过dyn1文件中(即Gamma点)的振动频率得到的**
**thermodynamics_from_dyns.csv 代表自由能是通过dyns文件中每个q点的振动频率得到的**
**thermodynamics_from_phtdos.csv 代表自由能是通过xxx_phono.dos声子态密度文件中对总态密度积分得到的。**


下面给出推到吉布斯自由能的推到过程。（参考链接：https://mp.weixin.qq.com/s?__biz=MzkyODI1MDkzMQ==&mid=2247483745&idx=1&sn=6d778b432d19423188077b703dd2dca1&chksm=c21aec3df56d652bc30408adaed234b6d1cc0bc7c23bfc816251ab5775cbbb999b6c8ab76cdb&scene=21#wechat_redirect）

一个材料体系，可以将其看作是由声子和电子组成。材料体系的总能当然是电子总能与声子总能之和。

电子总能包括电子间的相互作用能、原子核对电子的作用、范德华相互作用（如果考虑的话）等（这些都包含在第一性原理计算软件输出的能量中，以VASP为例，这一项就是使用“grep‘energy without entropy’  OUTCAR”这个命令搜索到的能量）以及电子熵对能量的贡献（这部分很小，较高温度、精确的计算中需要考虑，这无法直接从第一性原理计算软件中输出，而是需要额外的数据后处理）。

声子总能包括声子的内能以及声子的振动熵的贡献。而ZPE是声子总能的一部分。需要特别注意的是，第一性原理的计算软件通常只能够计算出声子的本征振动模（振动频率w），而声子的ZPE以及吉布斯自由能需要结合统计力学的知识处理得到。

既然ZPE是声子总能的一部分，那我们就从声子总能入手，看看ZPE到底是声子总能的那一部分？

对原子间的相互作用势能取简谐近似，并引入简正坐标消除势能交叉项，从而得到3N（N是体系中原子的个数）个谐振子的振动方程（详见：黄昆先生-固体物理学（1988年10月第一版p78-81））。根据谐振子的振动方程，求得谐振子的能量本征值（详见：钱伯初先生-量子力学（2006年1月第1版）p52）：

$$\varepsilon_{n_{\omega_{i}}} = (n_{\omega_i}+\frac{1}{2})\hbar\omega$$
$$i=1,2,3,...,3N, n_{\omega_i}=0,1,2,...,\infty$$

i等于3N，表示3N个谐振子。这3N个谐振子对应3N个声子，声子的频率分别为 W1....W3N（包含N个原子的体系共包含3N个声子）。将每个声子所有本征态看作一个独立的热力学系统，则这个声子本征态服从玻尔兹曼统计（下面涉及到的热力学公式详见：汪志诚先生-热力学.统计物理（第四版）第六、七章）。

玻尔兹曼统计的配分函数$Z_1$:

$$Z_1 = \sum_n e^{-\beta\varepsilon_n}$$
$$\beta=1/k_BT$$




###  <span style="color:yellow"> 电子态密度计算

####  <span style="color:green">**先进行非自洽计算获得dos数据, k点要比自洽密一些**</span>
```shell
scf -m mode=nscf execmd='mpirun -np 48' npool=4 queue=local kpoints_dense='16 16 16' 
```

####  <span style="color:green">**处理eletdos和elepdos数据获得可以origin绘图的数据**</span>
```shell
eletron -m mode=eledosdata execmd='mpirun -np 1' npool=1 queue=local kpoints_dense='16 16 16' 
```


算出来的*.tdos文件就是可以画图的总电子态密度数据了

第一列是能量值, 单位是eV

第二列是dos(E),根据$DOS(E) dE$定义, 表示在E和E+dE的能量范围内的能级数，单位是：**states/eV/(Unit Cell)**

注意：lambda.out中N(Ef)的单位是：**states/spin/Ry/(Unit Cell)**。如果想验证的话，可以从elph.inp_lambda.`number`文件获得相关信息：
```shell
     0.000000      0.000000      0.000000    10    60                            
     w1 w2 .... wn # n个振动模式
     Gaussian Broadening:   0.005 Ry, ngauss=   0     
     DOS = 11.556540 states/spin/Ry/Unit Cell at Ef= 16.274226 eV # 这里的DOS就是lambda.out中的 N(Ef)
     lambda(    1)=  0.0000   gamma=    0.01 GHz     
     lambda(    2)=  0.0000   gamma=    0.01 GHz
     ... ... ...
     ... ... ...
```

第三列是Int dos(E), 根据$\int_{E_0}^{E_1} DOS(E) dE$定义, 表示从最低能级到当前能级总的态的数量。

参考资料：http://hawk.fisica.uminho.pt/ricardo-ribeiro/QEnotes.html
```shell 
  E (eV)   dos(E)     Int dos(E) EFermi =   16.266 eV                                  
 -10.000  0.0000E+00  0.8000E+01
  -9.990  0.0000E+00  0.8000E+01
  -9.980  0.0000E+00  0.8000E+01
  ...     ...         ...
  ...     ...         ...
   9.980  0.2352E+01  0.3439E+02
   9.990  0.2339E+01  0.3442E+02
  10.000  0.2322E+01  0.3444E+02
```

###  <span style="color:yellow"> 电子能带结构计算

做完scf.in自洽计算，然后用pw.x做eleband.in的非自洽能带计算, 然后用band.x处理出能带的数据
```shell
# 很多人以为做qe的能带计算需要先做 calculation=’nscf’ 然后再做 calculation=’bands’
# 这其实是错误的, qe做能带计算只需要在 calculation=’scf’ 后做 calculation=’bands’ 即可.
The difference between a calculation=’bands’ and calculation=’nscf’ is that the former uses exclusively
the k-point provided while the latter might add points to respect crystal symmetry.
```

####  <span style="color:green">**获得eleband数据**</span>
```shell
eletron -m mode=eleband execmd='mpirun -np 1' npool=1 queue=local kinserted=200 nbnd=500
```

####  <span style="color:green">**处理eleband数据获得可以origin绘图的能带数据和投影能带数据**</span>
```shell
eletron -m mode=elebanddata execmd='mpirun -np 1' queue=local
```

#### <span style="color:green">**处理eleband的投影轨道数据获得可以获得origin绘制投影轨道的能带**</span>
```shell
eletron -m mode=elebandprojdata execmd='mpirun -np 1' queue=local
```

###  <span style="color:yellow"> 同时进行 电子能带结构计算 电子态密度计算 并且 处理好数据

<span style="color:green"> **(qe计算出来的电子态密度的横坐标能量的单位是Ry, 而vasp计算出来的电子态密度的横坐标的单位是eV)**

**1 Ry = 13.605693122994**

**1eV = 0.0734986443513 Ry**

```shell
eletron -m mode=eleproperties execmd='mpirun -np 48' npool=4 queue=local kinserted=200 nbnd=500  kpoints_dense='8 8 8' 
```
**不必处理数据，数据是已经都处理好的了，可以直接绘图**


###  <span style="color:yellow"> 计算超导

<span style="color:green"> **在使用lambda.x处理数据时，一定要注意3件事情**

<span style="color:green"> 1. 如果出现下面的错误，那么就需要修改lambda.f90的代码， nmodex改大一点
```shell
    wrong # or too many modes
```

```shell
#position: qe-7.1/PHonon/PH/lambda.f90
program elph
 
  ! read files 'filelph' produced by phonon (one for each q-point)
  ! sum over q-points to produce the electron-phonon coefficients:
  ! lambda (the one of BCS superconductivity) and alpha^2*F(omega)
  ! T_c using Allen-Dynes formula
 
  implicit none
  #integer, parameter:: npk=200, nsigx=50, nmodex=100, nex=200
  integer, parameter:: npk=200, nsigx=50, nmodex=600, nex=200
```


<span style="color:green"> 2. **如果elph_dir目录中有一个文件elph.inp_lambda.`number`文件中的频率是负值，那么就会导致lambda.out中无法处理处lambda和wlog.**

<span style="color:green"> 3. **如果分q点计算时，每一个q点目录中的scffit.in和scf.in用的都是相同的degauss值并且split_ph.in中的el_ph_sigma和el_ph_nsigma都相同，那么在每一个q点的elph.inp_lambda.`number`文件中都会有相同个数的DOS值并且DOS的值也都相同，这一点(即：DOS的值也都相同）至关重要.** 不然的话就会在lambda.out文件中报出一下错误：

```shell
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine lambda (4):
     inconsistent DOS(Ef) read
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     stopping ...
```

下面展示每个q点的elph_dir文件的elph.inp_lambda.`number`文件中的DOS at Ef的值：

```shell
# 例如使用`grep DOS elph.inp_lambda.1`，你会得到10个展宽对应的10个DOS值
grep DOS elph.inp_lambda.1 # 有多少个elph.inp_lambda.*文件被dyn0文件中的不可约q点数决定
     # 每个不可约q点有多少个DOS值被el_ph_nsigma决定，每个展宽的大小由l_ph_sigma和el_ph_nsigma共同决定
     DOS = 10.687167 states/spin/Ry/Unit Cell at Ef= 16.291748 eV #（DOS 和 Ef 值的大小被 scffit.in 和 scf.in 中的 degauss 影响）
     DOS = 10.436377 states/spin/Ry/Unit Cell at Ef= 16.288660 eV
     DOS = 10.257347 states/spin/Ry/Unit Cell at Ef= 16.287513 eV
     DOS = 10.299636 states/spin/Ry/Unit Cell at Ef= 16.287628 eV
     DOS = 10.470458 states/spin/Ry/Unit Cell at Ef= 16.289191 eV
     DOS = 10.651954 states/spin/Ry/Unit Cell at Ef= 16.291880 eV
     DOS = 10.824777 states/spin/Ry/Unit Cell at Ef= 16.295156 eV
     DOS = 11.028643 states/spin/Ry/Unit Cell at Ef= 16.298284 eV
     DOS = 11.315351 states/spin/Ry/Unit Cell at Ef= 16.300248 eV
     DOS = 11.712844 states/spin/Ry/Unit Cell at Ef= 16.299994 eV

# 例如使用`grep DOS elph.inp_lambda.2`，你会得到10个展宽对应的10个DOS值
grep DOS elph.inp_lambda.2
     DOS = 10.687167 states/spin/Ry/Unit Cell at Ef= 16.291748 eV
     DOS = 10.436377 states/spin/Ry/Unit Cell at Ef= 16.288660 eV
     DOS = 10.257347 states/spin/Ry/Unit Cell at Ef= 16.287513 eV
     DOS = 10.299636 states/spin/Ry/Unit Cell at Ef= 16.287628 eV
     DOS = 10.470458 states/spin/Ry/Unit Cell at Ef= 16.289191 eV
     DOS = 10.651954 states/spin/Ry/Unit Cell at Ef= 16.291880 eV
     DOS = 10.824777 states/spin/Ry/Unit Cell at Ef= 16.295156 eV
     DOS = 11.028643 states/spin/Ry/Unit Cell at Ef= 16.298284 eV
     DOS = 11.315351 states/spin/Ry/Unit Cell at Ef= 16.300248 eV
     DOS = 11.712844 states/spin/Ry/Unit Cell at Ef= 16.299994 eV
......
......
......
grep DOS elph.inp_lambda.10
     DOS = 10.687167 states/spin/Ry/Unit Cell at Ef= 16.291748 eV
     DOS = 10.436377 states/spin/Ry/Unit Cell at Ef= 16.288660 eV
     DOS = 10.257347 states/spin/Ry/Unit Cell at Ef= 16.287513 eV
     DOS = 10.299636 states/spin/Ry/Unit Cell at Ef= 16.287628 eV
     DOS = 10.470458 states/spin/Ry/Unit Cell at Ef= 16.289191 eV
     DOS = 10.651954 states/spin/Ry/Unit Cell at Ef= 16.291880 eV
     DOS = 10.824777 states/spin/Ry/Unit Cell at Ef= 16.295156 eV
     DOS = 11.028643 states/spin/Ry/Unit Cell at Ef= 16.298284 eV
     DOS = 11.315351 states/spin/Ry/Unit Cell at Ef= 16.300248 eV
     DOS = 11.712844 states/spin/Ry/Unit Cell at Ef= 16.299994 eV
```

**这两个文件中每一个展宽el_ph_sigma值对应的DOS和费米能级都一模一样。如果不一样的话，可能是你在算某一个q点时scffit.in和scf.in用的都是相同的degauss值与其它的不同导致。**

####  <span style="color:green"> 3. 解读声子态密度计算后得到的a2F.dos*文件。我们以H3S的a2F.dos1为例说明文件的含义

**H3S的原胞内有4个原子，所以其声子谱总共有12个振动模式。在a2F.dos1文件中除去第一列和第二列, 剩下的列数应该是12列，总共有14列。不信的话，你可以自己数一数。**

**所以：如果你想看每种元素对谱函数a2F的贡献，可以按照scf.in中元素顺序，依次求和，就可以获得每种元素的a2F。**

```shell
 
 # Eliashberg function a2F (per both spin)
 #  frequencies in Rydberg  
 # DOS normalized to E in Rydberg: a2F_total, a2F(mode) 

col1:频率      col2:总a2f    col3:H1       col4:H1        col5:H1       col6:H2        col7:H2        col8:H2       col9:H3       col10:H3       col11:H3       col12:S1       col13:S1       col14:S1
# 单位:里德伯Ry      无单位(后面列的同样没有单位)   
0.159727E-04    0.204760E-07    0.109697E-07    0.335915E-08    0.614710E-08    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00
0.479190E-04    0.552867E-06    0.296183E-06    0.907015E-07    0.165982E-06    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00
0.798652E-04    0.255958E-05    0.137122E-05    0.419918E-06    0.768447E-06    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00
0.111811E-03    0.702351E-05    0.376262E-05    0.115226E-05    0.210863E-05    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00
0.143758E-03    0.149275E-04    0.799694E-05    0.244898E-05    0.448162E-05    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00    0.000000E+00
......
......

```

####  <span style="color:green"> 4. 解读alpha2F.dat 

**还是以H3S为例说明。它是在计算lambda.in这一步时被处理出来的一个文件，如果你的动力学矩阵有虚频，那么这个文件的内容将全部是NaN.**

w_alpha2f_lambda.csv 可以用来绘制alpha-lambda的函数图像，它是在使用sc模块计算Tc时通过处理alpha2F.dat获得的。因为alpha2F.dat有多个展宽对应的a2f，所以需要根据实际情况选择一个收敛的展宽值。

```shell
# col1: 频率，单位时THz
# 之后的列都是各个展宽对于的总谱函数a2f. 展宽数以及展宽宽度取决于你的ph.in文件中参数设置
 E(THz)     0.005     0.010     0.015     0.020     0.025     0.030     0.035     0.040     0.045     0.050                                                                          
  0.0000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000  -0.00000
  0.3166  -0.00001  -0.00001  -0.00002  -0.00002  -0.00002  -0.00002  -0.00002  -0.00002  -0.00002  -0.00002
  0.6332  -0.00013  -0.00018  -0.00021  -0.00023  -0.00023  -0.00024  -0.00024  -0.00024  -0.00024  -0.00023
  0.9497  -0.00056  -0.00070  -0.00079  -0.00083  -0.00085  -0.00086  -0.00086  -0.00085  -0.00085  -0.00084
  1.2663   0.00020   0.00101   0.00135   0.00152   0.00162   0.00169   0.00172   0.00174   0.00174   0.00173
  1.5829   0.00584   0.00864   0.01015   0.01087   0.01129   0.01151   0.01158   0.01156   0.01148   0.01138
  1.8995   0.01644   0.01555   0.01631   0.01667   0.01675   0.01667   0.01648   0.01625   0.01600   0.01576
  2.2161   0.02681   0.01982   0.01865   0.01814   0.01745   0.01671   0.01602   0.01541   0.01490   0.01448
  2.5327   0.03587   0.03113   0.03035   0.03002   0.02910   0.02796   0.02690   0.02598   0.02521   0.02458
  2.8492   0.03952   0.03993   0.03987   0.03954   0.03866   0.03768   0.03682   0.03612   0.03555   0.03515
```


####  <span style="color:green"> 4. 解读ALPHA2F.OUT 

**ALPHA2F.OUT的第一列的单位是hartree，所以这里设计单位转化。有两种获得方法：**

**第1种：通过a2F.dos*获得：**
```shell
# a2F.dos5里面频率的单位是Ry, 转化为hartree需要除以2
file=a2F.dos5
sed '1,5d' $file | sed '/lambda/d' | awk {print $1/2, $2} > ALPHA2F.OUT
```

**第2种：通过alpha2f.dat获得：**
```shell
# alpha2f.dat里面频率的单位是Thz, 转化为hartree需要除以6579.684
# $6代表选择第6列也就是第5个展宽（因为第一列是频率）
sed '1,1d' alpha2f.dat | awk '{print $1/6579.684, $6}' ALPHA2F.OUT
```

```shell
# 第一列代表频率，单位是hartree
# 第二列代表总谱函数a2f,
9.306e-06 0.363609E-06
2.82051e-05 0.983660E-05
4.71043e-05 0.455576E-04
6.60035e-05 0.125031E-03
8.49025e-05 0.265761E-03
0.000103801 0.485253E-03
0.000122701 0.801010E-03
0.0001416 0.123654E-02
0.000160499 0.189306E-02
0.000179399 0.302340E-02
```

**处理电荷屏蔽常数为0.1和0.13,得到 $\lambda$ 和 $\omega_{log}$ 并且输出一个文件**
```shell
# 根据真实费米面选择
sc -m mode=Tc execmd='mpirun -np 48' npool=1 queue=local temperature_steps=100 a2fdos=True alpha2fdat=False broaden=0.5 smearing_method=1 nef=xxxxx
```

```shell
# 自动选择gaussid以及gauss
sc -m mode=Tc execmd='' npool=1 queue=local temperature_steps=100 a2fdos=True alpha2fdat=False broaden=0.5 smearing_method=1
```

```shell
# 当然如果你不满意程序给你选择的gaussid以及gauss, 你也可以自己指定
qe_main.py -i relax.out -j bash sc -m mode=Tc execmd='mpirun -np 48' npool=1 queue=local temperature_steps=100 a2fdos=True alpha2fdat=False broaden=0.5 smearing_method=1 gaussid=10 gauss=0.05
```


```shell
# 指定最高频率
top_freq=80

# 当然了，a2fdos， alpha2fdat 两个参数不设置也可以，
# 默认使用的是从alpha2fdat中读取数据，
# 因为很多时候，你不会计算phonodos，所以你也没有a2F.dos*这些文件。

# 使用eliashberg方法超导转变温度, 指定读取a2F.dos*文件
a2fdos=True alpha2fdat=False

# 使用eliashberg方法超导转变温度, 
# 指定读取alpha2F.dat文件中使用哪一列的degauss对应的alpha2F数值。
# 使用alpha2fdat来指定
# 这个方法生成ALPHA2F.OUT可能有问题导致 ELIASHBERG_GAP_T.OUT 中出现NAN。 所以更推荐a2fdos=True那种处理方法。
a2fdos=False alpha2fdat=True
```

**(这个方法生成ALPHA2F.OUT可能有问题导致 ELIASHBERG_GAP_T.OUT 中出现NAN。所以更推荐a2fdos=True那种处理方法。)**

用到的公式：

$$\lambda = 2 \int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega$$

$$\omega_{log} = exp[\frac{2}{\lambda}\int_{0}^{\infty} \frac{d\omega}{\omega} \alpha^2 F(\omega) ln(\omega)]$$

####  <span style="color:green"> **千万注意：alpha2F.dat中频率的单位是THz, 但是freq.gp文件中的频率的单位是cm-1，如果想在一张图中把两个数据放在一起对比需要一个单位转化**
$$c=299792458  m/s$$
$$\lambda^{-1} = \frac{\nu}{c} = \frac{1Thz}{299792458  m/s} =  \frac{10^{12}Hz}{299792458\times 10^{2} cm \cdot Hz} = 33.3564095198152 cm^{-1} $$
即：
$$1Thz \Leftrightarrow 33.3564 cm^{-1}$$

w_alpha2f_lambda.csv 中的第一列是频率，其单位是经过转化的$cm^{-1}$



####  <span style="color:green"> **Eliashberg方法计算Tc需要INPUT文件中只需设置两个参数**
1. 前者是screen_constant，一般取0.10~0.13；
2. 后者是temperature_steps，表示对ntemp个温度点自洽求解Eliashberg方程，得到带隙Δ关于温度T的曲线(该程序首先处理McMillan方程，得到超导临界温度tc作为参考值，然后在温度区间[tc/6, tc*3]中线性插入ntemp个温度点。ntemp一般取40~100即可，也可以更大，建议根据体系差异灵活调控)。

####  <span style="color:green"> **Eliashberg方法计算得到的ELIASHBERG_GAP_T.OUT文件就是能隙方程关于温度的关系，第一列是温度，第二列是能隙（能隙的单位是hartree, 所以如果你想要转化为meV，就需要: 第二列 * 27211.38602)**
1 Hartree = 27.21138602 eV 

1 eV = 0.03674932 Hartree

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

###  <span style="color:yellow"> 计算动力学矩阵元，获得高对称路径下的声子振动频率和声子线宽

<span style="color:green"> **这一步会覆盖上一步计算声子态密度时产生的文件gam.lines 和 xxx.freq.qp 和 xxx.freq 三个文件。所以一定要确保计算完声子的态密度之后将声子态密度数据处理出来。**

```shell
phono -m mode=matdyn qpoints='6 6 6' execmd='mpirun -np 48' npool=1 queue=local qinserted=50
```

####  <span style="color:green"> 处理数据获得声子谱
<span style="color:green"> **注意！！！ 计算完高对称路径下的声子振动频率后，一定要记得处理数据获得声子谱。输出一个文件 qp_freq_width.csv 可以用来绘制声子谱以及每个振动频率对应的线宽**

####  <span style="color:green"> 这里是处理声子谱的命令

```shell
phono -m mode=phonobanddata execmd='mpirun -np 48'
```
####  <span style="color:green"> **这里是获得高对称点投影路径的命令， 可以用来在origin中绘制高对称点的垂直参考线**
```shell
phono -m mode=hspp execmd='mpirun -np 1'
```

**因为计算声子谱设置的q点和计算声子态密度设置的q点不一样，但是它们两个计算完都会重新覆盖，`体系名称`.freq文件, `体系名称`.freq.gq文件这三个文件。所以请务必确保：无论是计算完声子态密度还是声子谱，都一定要处理数据，然后再进行下一步计算。**

**不然，可能出现声子的q点数目就和gam.lines文件中q点数目对不上的情况。**

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




###  <span style="color:yellow">  自动prepare计算
**只需要指定qpoints即可，kpoints_dense和kpoints_sparse会自动通过qpoints*4， qpoints*2来获得**
```shell
prepare -m mode=prepareall electron_maxstep=1000 core=4 npool=1 queue=local qpoints='4 4 5'
```
```shell
prepare -m mode=preparescf electron_maxstep=1000 core=4 npool=1 queue=local  qpoints='4 4 5'
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

#### 最繁参数
```shell
relax -m mode=rvf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 三次结构弛豫  </span>

#### 最繁参数
```shell
relax -m mode=rv3 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 单次结构弛豫 </span>

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

#### 最繁参数
```shell
batchrelax -m mode=rvf core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 批量结构弛豫 同时 每个结构做三次结构弛豫  </span>

#### 最繁参数
```shell
batchrelax -m mode=rv3 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

### <span style="color:yellow"> 批量结构弛豫 同时 每个结构做三次结构弛豫  </span>

#### 最繁参数
```shell
batchrelax -m mode=rv1 core=28 ediff=1e-8 ediffg=-0.001 ismear=1 kspacing=0.18 encut=800
```

###  <span style="color:yellow"> 有限位移法计算声子谱  </span>

#### 最繁参数
```shell
# 这里有两种设置KPOINTS的方法：
# 方法一：根据kdensity='36 36 36'进行设置，命令如下：
phono -m supercell='2 2 2' kdensity='36 36 36' mode=disp execmd='mpirun -np 48' ismear=1 encut=800 ediff=1E-08 ediffg=-0.001 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'
# 方法二：根据kspacing=0.18进行设置，命令如下：
phono -m supercell='2 2 2' kspacing=0.18 mode=disp execmd='mpirun -np 48' ismear=1 encut=800 ediff=1E-08 ediffg=-0.001 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'
```

###  <span style="color:yellow"> 密度泛函微扰DFPT法计算声子谱  </span>

#### 最繁参数
```shell
# 这里有两种设置KPOINTS的方法：
# 方法一：根据kdensity='36 36 36'进行设置，命令如下：
phono -m supercell='2 2 2' kdensity='36 36 36' mode=dfpt execmd='mpirun -np 48' ismear=1 encut=800 ediff=1e-08 ediffg=-0.001 lreal=.FALSE. queue=lhy
# 方法二：根据kspacing=0.18进行设置，命令如下：
phono -m supercell='2 2 2' kspacing=0.18 mode=dfpt execmd='mpirun -np 48' ismear=1 encut=800 ediff=1e-08 ediffg=-0.001 lreal=.FALSE. queue=lhy
```

###  <span style="color:yellow"> 有限位移法计算声子谱——数据处理band  </span>

#### 最繁参数
```shell
data -m mode=dispband supercell='2 2 2' spectrum=True 
```

###  <span style="color:yellow"> 有限位移法计算声子谱——数据处理phdos  </span>

#### 最繁参数
```shell
data -m mode=dispphdos supercell='2 2 2' spectrum=True 
```

###  <span style="color:yellow"> 密度泛函微扰DFPT法计算声子谱——数据处理band  </span>

#### 最繁参数
```shell
data -m mode=dfptband supercell='2 2 2' spectrum=True
```

###  <span style="color:yellow"> 密度泛函微扰计算声子谱——数据处理phdos  </span>

#### 最繁参数
```shell
data -m mode=dfptphdos supercell='2 2 2' spectrum=True 
```


###  <span style="color:yellow"> 电子性质计算  </span>

如果通过-w指定工作路径，那么就在指定路径下进行电子性质计算。

如果没有通过-w指定路径，那么就在当前路径下进行电子性质计算。

当前可以计算的电子性质包括：

    1.自洽: scf 
        在指定的目录或默认的目录下创建scf目录，例如 -w 指定了eletron，那么实际自洽计算的目录为eletron/scf
        1.1 在计算自洽时同时打开了LCHARG =.TRUE.和LAECHG =.TRUE.，所以会同时计算bader电荷转移
        1.2 在计算自洽时同时打开了LELF =.TRUE.，所以会同时计算ELF电子局域函数
    2.电子态密度: eledos
        2.1 在指定的目录或默认的目录下创建eledos目录，例如 -w 指定了eletron，那么实际自洽计算的目录为eletron/eledos。计算eledos时需要将scf的CHGCAR拷贝过来。
        2.2 电子态密度计算时，其kspacing需要是scf计算的2倍
    3.能带: eband
        在指定的目录或默认的目录下创建eband目录，例如 -w 指定了eletron，那么实际自洽计算的目录为eletron/eband。计算eledos时需要将scf的CHGCAR拷贝过来。
    4. 晶体轨道重叠布居COOP和COHP
        参考网页: 
            https://blog.shishiruqi.com//2019/04/22/cohp/
            https://zhuanlan.zhihu.com/p/470592188
        参考书：
            Roald Hoffman-Solids and Surfaces_ A Chemist’s View of Bonding in Extended Structures (1988)
        参考文献：
            J. Phys. Chem. 1993, 97, 8617-8624
            J. Phys. Chem. A 2011, 115, 5461-5466
            J. Comput. Chem. 2016, 37, 1030-1035

**电子态密度计算时，INCAR中有几个参数需要格外注意，这里我将这些参数给出，请你使用该脚本产生eledos计算的INCAR后稍作检查(原因是我在设置了EMIN和EMAX后算出来的dos的能量范围没有正值。)**
```shell
ISTART = 1    # 读取WAVECAR，如果WAVECAR不存在，就自己构造
ICHARG = 11   # 读取CHGCAR， 进行非自洽计算
#EMIN = -10   # 该脚本中默认将其注释了，EMIN指定了DOS评估的能量范围的下限。
#EMAX =  10   # 该脚本中默认将其注释了，EMAX指定了DOS评估的能量范围的上限。
ISMEAR  = -5  # For DOS
LORBIT = 11   # 输出分波态密度信息
```

#### 计算命令介绍

```shell
# 只计算scf,  工作路径-w下一定包含eband, scf, eledos三个文件，即使单独算其中一个，也要在eletron目录下，不要进入scf目录计算。
eletron -m mode=scf encut=800 ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


# 只计算能带,  工作路径-w下一定包含eband, scf, eledos三个文件，即使单独算其中一个，也要在eletron目录下，不要进入eband目录计算。在准备输入阶段，会先检查eletron目录中有没有scf/CHGCAR，如果没有就退出。
eletron -m mode=eband  encut=800 ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


# 只计算电子态密度,  工作路径-w下一定包含eband, scf, eledos三个文件，即使单独算其中一个，也要在eletron目录下，不要进入eledos目录计算。在准备输入阶段，会先检查eletron目录中有没有scf/CHGCAR，如果没有就退出。
# 前面提到自洽的kspacing是电子态密度的一半，不需要你自己指定，只需要保持你设置的kspacing与scf的kspacing保持一致即可，程序自己会加倍kspacing。
eletron -m mode=eledos  encut=800 ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


# 同时计算scf和eband
eletron -m mode='scf eband'  encut=800 ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


# 同时计算scf和eledos。# 前面提到自洽的kspacing是电子态密度的一半，不需要你自己指定，只需要保持你设置的kspacing与scf的kspacing保持一致即可，程序自己会加倍kspacing。
eletron -m mode='scf eledos' encut=800  ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


# 同时计算scf和eledos和eband。# 前面提到自洽的kspacing是电子态密度的一半，不需要你自己指定，只需要保持你设置的kspacing与scf的kspacing保持一致即可，程序自己会加倍kspacing。
eletron -m mode='scf eledos eband' encut=800  ediff=1e-8 ediffg=-0.001 ismear=0 sigma=0.01 kspacing=0.18 lreal=.FALSE. queue=lhy execmd='mpirun -np 48'


#获得投影到一维路径的高对称路径点
vasp_main.py -i CONTCAR -j bash data -m mode=hspp execmd='mpirun -np 48'
```



#### 常见错误分析
```shell
错误: charge density could not be read from file CHGCAR for ICHARG>10

分析: 可能原因是scf的NGX, NGY, NGZ与eband的不同, 但是归根结底还是可能设置的kpoints出了问题。

解决: 可以通过手动指定 NGX,NGY,NGZ in eledos/INCAR 与 scf/OUTCAR 的保持一致。
```
```shell
错误: Your FFT grids (NGX, NGY, NGZ) are not sufficient for an accurate calculation. Thus, the results might be wrong. 

分析: 可能原因是scf的NGX, NGY, NGZ与eband的不同,有可能是因为PREC设置的精度在scf和eband中不同导致。

解决: 可以通过手动指定 eband/INCAR 的 PREC 与 scf/INCAR保持一致
```
```shell
警告: dimensions on CHGCAR file are different

错误: charge density could not be read from file CHGCAR for ICHARG>10

分析：出现这个这个错误的可能原因有：
1. 自洽的结构和能带、dos的结构不一致
2. 自洽的INCAR的PREC和ENCUT 与 能带、dos的INCAR的PREC和ENCUT不一致
```

```shell
最近在计算有限位移的声子谱和自洽测encut收敛性时遇到一个错误

警告: PSMAXN for non-local potential too small

错误: REAL_OPTLAY: internal error (1)

分析: 这可能是由于LREAL=.TRUE. 或者 LREAL=Auto导致.两个办法解决：
    1. 降低ENCUT，可以保持LREAL参数不变让vasp收敛
    2. ENCUT不变，关闭LREAL

分析：出现这个这个错误的可能原因有：
1. 自洽的结构和能带、dos的结构不一致
2. 自洽的INCAR的PREC和ENCUT 与 能带、dos的INCAR的PREC和ENCUT不一致
```





# <div align="center"> <span style="color:red"> mytoolkit篇 </span> </div>

## 格式转化
```shell
tool_main.py -i 输入文件名称 -w ./ convert -m dst_format=输出文件名称
输出文件名称必须包含后缀为.vasp, .cif, .struct的内容
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