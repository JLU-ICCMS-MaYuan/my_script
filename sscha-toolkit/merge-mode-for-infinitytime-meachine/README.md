#  <span style="color:orange"> merge-mode-for-infinitytime-meachine

这种模式是针对可以无限占用节点的机器开发出来的计算模式, 也是sscha官方非常推崇的方式.

##  <span style="color:red"> 0. 安装sscha

由于sscha提交任务的模式仅支持通过ssh连接到远程, 将任务分发到远程机器, 然后打包回收到本地做结构弛豫. 所以无法在主节点上运行该任务.
为了在主节点运行该任务, 直接使用slurm提交任务, 需要修改sscha的源码.

将该仓库中的Cluster.py整个文件拷贝到sscha环境中相应的源码位置即可. 如何找到相应的sscha环境中的源码呢?
一般来说, 先激活conda 环境`conda activate sscha`. 然后进入sscha环境. `cd  ~/miniconda/env/sscha`这是我的路径. 然后进入`cd  lib/python3.11/site-packages/sscha`就找到Cluster.py了.

另一个非常重要的问题：就是Julia中调用python相关库的安装问题。sscha中，提到了
```shell
You must ensure Julia’s dependencies are correctly set up to activate Julia’s speedup on the SSCHA minimization. To do this, run the following line:

python -c 'import julia; julia.install()'
Note: this command may fail if you are using micromamba. To solve the issue, you need to manually specify the binary location of micromamba to julia:

export CONDA_JL_CONDA_EXE=$HOME/.local/bin/micromamba
Replacing $HOME/.local/bin/micromamba with the path to the micromamba binary if you changed the default. To make it work after the next login, add the environment variable to the init script

echo "export CONDA_JL_CONDA_EXE=$HOME/.local/bin/micromamba" >> $HOME/.bashrc
To configure Julia PyCall to work with anaconda (or micromamba), open a Julia shell, typing Julia. Enter in the package manager by typing ]. You should see your prompt turning into a pkg>. Then build the conda extension and compile PyCall.

build Conda
add PyCall
```

你在激活了Julia之后，用`build Conda`和`add PyCall`就可以安装相关的库。但是在实际操作中，我发现这两个命令都无法成功安装。下面的命令是可以使用的：
```julia
using Pkg # 加载 Pkg 模块
Pkg.build("Conda")
Pkg.add("PyCall") 
```


## <span style="color:red"> 1. 准备初始动力学矩阵

做一个粗q网格的声子计算(比如2x2x2),获得动力学矩阵. 在0.generate-init-dyn目录下进行.
获得粗网格的动力学矩阵后，执行0.generate-init-dyn目录下面的copyDYNtoDYN_POP.sh
命令如下：
```shell
. copyDYNtoDYN_POP.sh 4 MgB2
```

## <span style="color:red"> 2. 将0.generate-init-dyn中准备的初始动力学矩阵拷贝到1.relax目录下, 并修改relax_withQE.py目录下的relax.py文件中的以下参数:

需要修改relax_withQE.py中的参数. (relax_withQE.py是所有的自洽计算用QE计算, relax_withMTP.py是所有的自洽计算都用mtp计算) 


```python

# 0.frequently changed parameters
BATCH_SIZE = 20     # 总共提多少个任务
JOB_NUMBER = 5      # 每个脚本中包含多少个自洽计算
ONSET_DYN_POP_IDX = 0 # 你准备读入的动力学矩阵的名称, 他在relax函数的start_pop参数也用到了, 但是为dyn_pop_idx+1: 这意味着要把新生成的动力学矩阵编号记为dyn_pop_idx+1
MAX_POPULATION = 10   # 算到第多少代终止
TARGET_PRESSURE = 200 # 指定压强 GPa


# 1.Prepare QE input parameters部分

atom_masses = {}
pseudo = {}
input_data = {pseudo_dir, ecutwfc, degauss, conv_thr, mixing_beta}
k_points = ()

# 2.Load the dynamical matrix
irr_idx # 不可约动力学矩阵的个数

# 3.Prepare random configurations
T0

# 4.Rrepare the cluster
mpi_cmd
AlreadyInCluster # AlreadyInCluster=True代表本地提交任务, 既可以把当前python脚本提交到节点上, 也可以在主节点运行当前脚本 AlreadyInCluster=False代表通过ssh连接远程机器提交任务
my_hpc.workdir # 大量运行qe自洽计算的路径
my_hpc.binary # pw.x的路径和运行方式
my_hpc.load_modules
my_hpc.n_cpu
my_hpc.n_nodes
my_hpc.n_pool
my_hpc.set_timeout(600) # 每个任务限制10min时长
my_hpc.time = "100:00:00" # 每个脚本限制100hours时长

```

## <span style="color:red"> 3. 迭代好动力矩阵后, 使用CalHess.py计算hessian矩阵
### <span style="color:yellow"> 3.1 修改计算脚本 CalHess.py, 有以下 10 个参数需要修改
```python
DATA_DIR="data_ensemble_manual" # 从该目录中获得动力学矩阵和能量受力, 有时候这个参数要修改为popj_N, j代表代数, N代表该代总计随机结构数(例如: pop12_2000代表从第12代的随机结构的能量受力搞出V3_Hessian.dyn)
N_RANDOM = 1000 # 动力学矩阵达到收敛的那一代使用的总随机结构的数量
DYN_PREFIX =  'dyn_pop11_' # 达到收敛的那一代动力学矩阵的前一代动力学矩阵
FINAL_DYN =   'dyn_pop12_' # 达到收敛的那一代动力学矩阵
SAVE_PREFIX = 'V3_Hessian.dyn' # 计算出的Hessian矩阵的结果
NQIRR = 4 # 总的动力学矩阵的数量
Tg = 0 # 温度
T =  0 # 温度
POPULATION = 12 # 达到收敛的那一代动力学矩阵的编号, 程序从DATA_DIR="data_ensemble_manual"中读取那一代动力学矩阵计算出的受力和能量
ens.load(DATA_DIR, population = POPULATION , N=N_RANDOM ) # 注意这里load函数加载的文件格式是.dat, 如果DATA_DIR目录中存放的文件格式是.npy, 
#ens.load_bin(DATA_DIR, population_id = POPULATION), 特别的load_bin函数没有 N=N_RANDOM 这个参数.
```

在提交完计算任务后，可以通过`python check_sscha_efv.py sscha.out`命令检查计算是否收敛。
如果做的是NPT系综，需要检查FC, gibbs, STRESS TENSOR是否收敛。即：FC是否收敛到1e-4数量级，gibbs是否收敛到1meV/atom，STRESS TENSOR是否收敛到指定压强附近0.1GPa.

如果做的是NVT习总，只需要同理需要关注上述三个物理量的，但是STRESS TENSOR只要稳定到程序自动计算的压强点即可。


### <span style="color:yellow"> 3.2 如果计算完成收敛，通过CalHess.py计算出HESS矩阵

将CalHess.py提交到节点运行, 可能需要你根据自己使用的集群修改提交作业的脚本
```shellV
sbatch SUB_HESS.sh
```
计算完成后会得到一系列的`V3_Hessian.dyn*`，需要检查其中是否有虚频，最好是没有虚频。如果有虚频，可能是计算有误，也可能是结构本身不稳定。

### <span style="color:yellow"> 3.2 将1_CalHess.py提交到节点运行, 可能需要你根据自己使用的集群修改提交作业的脚本
```shell
sbatch SUB_HESS.sh
```


## <span style="color:red"> 4. 读取V3_Hessian.dyn1（V3_Hessian.dyn*内部的所有结构信息都是相同的，挑选第一个即可）矩阵中的晶格参数和原子坐标, 在2.interpolation/1.sparse中计算粗的动力学矩阵, 在2.interpolation/2.fine里面计算细的动力学矩阵


### <span style="color:yellow"> 4.1 重新计算稀疏的、稠密的动力学矩阵

需要重新准备自洽和声子计算的输入文件

输入文件`scffit.in`和`scf.in`中的结构信息最好都从`V3_Hessian.dyn`中提取。

`V3_Hessian.dyn`中的结构文件的说明:
```shell
Dynamical matrix file
File generated with the CellConstructor by Lorenzo Monacelli
#  9.4567007000000007 是 celldm(1)
3 27 0     9.4567007000000007     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000
Basis vectors # 单位是alat, 对应的scffit.in scf.in中都要用alat
    1.0092654065722451    -0.0000000000000000     0.0000000000000000
   -0.5046327032861226     0.8740494812534840     0.0000000000000000
    0.0000000000000000     0.0000000000000000     0.6998060393036162
        1  'Ce '  127707.9215674129955005
        2  'Sc '  40974.8071860994023154
        3  'H '    918.6811103989390404
    # 单位也是alat
    1     1    -0.0000000000000000     0.0000000000000000    -0.0000000000000000
    2     2    -0.0000000000000001     0.5826996540344207     0.3499030198531194
    3     2     0.5046327032861220     0.2913498270172102     0.3499030198531194
    4     3     0.6227139914744079    -0.0000000000000002     0.5498098458271082
    5     3    -0.3113569957372042     0.5392861358209859     0.5498098458271082
    6     3     0.1932757075489183     0.3347633452306448     0.5498098458271082
    7     3     0.3113569957372040     0.5392861358209854     0.1499961938791310
    8     3     0.3865514150978366    -0.0000000000000002     0.1499961938791310
    9     3    -0.1932757075489188     0.3347633452306448     0.1499961938791310
   10     3     0.2303073226102307    -0.0000000000000000     0.3499030198531194
   11     3    -0.1151536613555787     0.1994519920201268     0.3499030198531194
   12     3     0.3894790419305433     0.6745974890315040     0.3499030198531194
   13     3     0.1151536613555785     0.1994519920201264     0.3499030198531194
   14     3     0.7789580839620139    -0.0000000000000001     0.3499030198531194
   15     3    -0.3894790419305441     0.6745974890315038     0.3499030198531194
   16     3     0.3865514150978366    -0.0000000000000002     0.5498098458271082
   17     3    -0.1932757075489188     0.3347633452306448     0.5498098458271082
   18     3     0.3113569957372040     0.5392861358209854     0.5498098458271082
   19     3     0.1932757075489183     0.3347633452306448     0.1499961938791310
   20     3     0.6227139914744079    -0.0000000000000002     0.1499961938791310
   21     3    -0.3113569957372042     0.5392861358209859     0.1499961938791310
   22     3     0.5046327032861222     0.4563658209178845     0.0000000000000000
   23     3     0.3617246604406449     0.2088418300668730     0.0000000000000000
   24     3     0.6475407461315990     0.2088418300668730     0.0000000000000000
   25     3     0.0000000000000000     0.4176836601337461    -0.0000000000000000
   26     3     0.1429080428454766     0.6652076509847578     0.0000000000000000

```

如果不愿意手动制作`scffit.in`, `scf.in`的输入文件，可以用`qedyn2struct.py`脚本制作。执行命令如下。注意：最好去读取`V3_Hessian.dyn*`中的结构.
```shell
# 将其转化为POSCAR
python qedyn2struct.py -i V3_Hessian.dyn1 -o vasp

# 将其转化为qe的scf.in输入文件
python qedyn2struct.py -i V3_Hessian.dyn1 -o qe -kd 12 12 12 -ks 6 6 6 -pp 赝势的目录 -pn Mg Mg_uspp.upf B B_uspp.upf
```



### <span style="color:yellow"> 4.2 回收`1.sparse`目录中粗q网格的动力学矩阵，`2.fine`目录中细q网格的动力学矩阵，`1.relax`目录中V3_Hessian动力学矩阵。
用到了`3.Inter`里面的`get_dyn`脚本.


如果报错，`D_S (l=2) for this symmetry operation is not orthogonal`，
这说明V3_Hessian.dyn可能弛豫出来的对称性微微扭曲，建议在那一代补充增加结构数，获得新的V3_Hessian.dyn矩阵。
```shell
task #         0
from d_matrix : error #         3
D_S (l=2) for this symmetry operation is not orthogonal
```

### <span style="color:yellow"> 4.2 回收`1.sparse`目录中稀疏q网格的动力学矩阵，`2.fine`目录中稠密q网格的动力学矩阵，`1.relax`目录中V3_Hessian动力学矩阵。
用到了`3.Inter`里面的`get_dyn`脚本.
```shell
. get_dyn.sh <sparse_qmesh> <sparse_dyns> <fine_qmesh> <fine_dyns> <prefix>
# 例如：
. get_dyn.sh 222 4 666 28 MgB2
```



## <span style="color:red"> 5. 基于2.interpolation/1.sparse中稀疏的动力学矩阵和2.interpolation/2.fine中稠密的动力学矩阵, 使用脚本2_Inter.py在2.interpolation/3.Inter中计算插值的动力学矩阵


修改`2_Inter.py`的内容，并通过`SUB_INTER.sh`提交任务。
```python
···
...
# 详细说明2_Inter.py需要修改的内容：
inter_dyn = CC.Phonons.Phonons()
dyn_sscha = CC.Phonons.Phonons("V3_Hessian.dyn",4) # 指定好V3_Hessian的名称和不可以q点数量。
dyn_coarse= CC.Phonons.Phonons("222.dyn",4) # 指定好粗动力学矩阵的名称和不可以q点数量。
dyn_fine  = CC.Phonons.Phonons("666.dyn",20)# 指定好细动力学矩阵的名称和不可以q点数量。
# 以上名称有赖于你自己的定义。前后统一即可。
inter_dyn = dyn_sscha.Interpolate(coarse_grid=[2,2,2], fine_grid=[6,6,6], support_dyn_coarse=dyn_coarse, support_dyn_fine=dyn_fine, symmetrize=True)
#dyn_222.SwapQPoints(dyn_sscha)
#dyn_222.save_qe("i")
inter_dyn.save_qe("inter_dyn_") # 插好值的动力学矩阵被命名为inter_dyn_1, inter_dyn_2, ......
···
···

```

```shell
sbatch SUB_INTER.sh
```

## <span style="color:red"> 6. 基于插值的动力学矩阵inter_dyn_*, 在3.Tc目录中计算电声耦合.
### <span style="color:yellow"> 6.1 首先，将`inter_dyn_*`, 降低拷贝为指定体系名称的动力学矩阵。只需要看看你之前的动力学矩阵的命名规则就行。
### <span style="color:yellow"> 6.2 其次在`3.Tc`目录中准备如下四个文件：`scffit.in`(自洽输入文件la2F=.true.), `scf.in`(自洽输入文件la2F=.false.), `s5_PhAssignQ.sh`(提交任务的脚本), `split_ph.in`(分q点计算的脚本)。特别注意：
   ```shell
   Electron-phonon coefficients for Nb4H14
    &inputph
    tr2_ph=1.0d-14,
    prefix='Nb4H14',
    fildvscf='dvscf',
    electron_phonon='interpolated',
    el_ph_sigma=0.005,
    el_ph_nsigma=10,
    alpha_mix(1)=0.3,
    amass(1)=92.90638 ,
    amass(2)=1.00794 ,
    outdir='./tmp',
    fildyn='Nb4H14.dyn',
    trans=.false., # 一定要false，代表从动力学矩阵后续算电声耦合计算
    ldisp=.false., # 一定要false，代表通过指定q点坐标的方式进行分q点计算
    nq1=6,nq2=6,nq3=6,
    start_q=1
    last_q=1
    /
    XQ1 XQ2 XQ3
   ```

```shell
   Error in routine initialize_grid_variables (1):
     problems reading u
```


### <span style="color:yellow"> 6.3 从`2.fine`目录中拷贝记录了不可约q点坐标的`*.dyn0`文件到`3.Tc`目录中。
   
### <span style="color:yellow"> 6.4 准备分q点目录，并且在各个分q点目录中tmp文件中的`dv文件`和`patterns文件`。
   
    **这里非常重要。因为通过`ldisp=.false.`和`XQ1 XQ2 XQ3`方式指定q点，所以每一个分q点的tmp/_ph0/*.phsave/目录中patterns.*.xml都要改为tmp/_ph0/*.phsave/patterns.1.xml**

    **<span style="color:lightblue"> 以上步骤都可以通过`pre.sh`准备好, 准备好之后通过`sub.sh`提交任务。**

5. 后面的步骤就和`QE`计算超导的步骤一样了。

## <span style="color:red"> 7. 报错集锦

### <span style="color:yellow"> 7.1 一个很棘手的问题：自洽不收敛。

当你的体系中有很多电子时候，qe自动生成的nbnd就不太够用了。你需要自己大概计算一下。

**怎么自己计算nbnd呢？**

绝缘体：nbnd = 价带数目 = 总电子数/2

金属：nbnd = 总电子数/2 * 1.2， 最小为4

例如：CeH9是一个金属，里面包含2个Ce和18个H. 每个Ce贡献12个电子，每个H贡献1个电子，总电子数为12*2+18=42. 所以nbnd = 42/2 * 1.2 = 25.2. 所以四舍五入 nbnd=25. 这与qe自动生成的nbnd=25是一致的.
```shell
# 这是CeH9计算自洽scf.out中的一部分输出，可以看到总共有25条带。25条带中有23条被占据。
          k = 0.0000 0.0000 0.0000 (  5347 PWs)   bands (ev):

   -17.5304 -16.8818  -4.0034  -0.0866  -0.0866   0.2847   0.2847   0.6015
     1.4679   4.3284   8.8980   9.6431   9.6431  10.0894  11.2221  11.2221
    12.0650  12.1232  13.4139  13.4139  13.6684  14.2188  14.2188  18.3537
    18.4474

     occupation numbers
     1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000
     1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000
     1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   0.0000
     0.0000
```

这是单胞CeH9的nbnd=25, 如果是2x2x2的CeH9，里面包含了160个原子，总共空带需要25*8=200条。但是你会发现此时自洽计算会不收敛。
你需要大概设置比200再多一点，至于多多少，需要自己试一试。激进一点，直接写300我觉得问题也不是很大。


### <span style="color:yellow"> 7.2 续算是结构优化失败。

很可能是因为随机结构太少，或者随机结构不合理，无法优化。可以尝试多产生一点随机结构。

```shell
ERROR WHILE UPDATING THE WEIGHTS

Error, one dynamical matrix does not satisfy the acoustic sum rule.
    If this problem arises on a sscha run,
    it may be due to a gradient that violates the sum rule.
    Please, be sure you are not using a custom gradient function.

DETAILS OF ERROR:
    Number of translatinal modes in the original dyn = 0
    Number of translational modes in the target dyn = 3
    (They should be both 3)

```