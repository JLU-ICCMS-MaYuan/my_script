# 流程大纲

## <span style="color:red"> 1. 根据第j代的动力学矩阵, 产生随机结构, 生成随机结构自洽计算的输入文件, 生成提交任务的脚本, 提交任务, 检查每个随机结构自洽计算的结果, 回收自洽计算的结果, 结构弛豫并生成第j+1代动力学矩阵

### <span style="color:yellow"> 1.1 dyn文件改名: 222q网格计算出来的4个dyn文件分别改名为
```shell
mv Ce1Sc2H24.dyn1 dyn_pop0_1
mv Ce1Sc2H24.dyn2 dyn_pop0_2
mv Ce1Sc2H24.dyn3 dyn_pop0_3
mv Ce1Sc2H24.dyn4 dyn_pop0_4
```
dyn_pop0_1 代表第0代的第1个动力学矩阵，第0代也就是初始代

### <span style="color:yellow"> 1.2 生成随机结构1_CreatEns.py
1_CreatEns.py是生成自洽文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R1_CreatEns.py, 通过运行R1_CreatEns.py生成随机结构。

R1_CreatEns.py将生成的随机结构相关的文件存放在data_ensemble_manual中，分别是

```
scf_population1_1.dat, ..., scf_population1_100.dat
scf_population1_1.dat 代表第1代随机结构的第1个随机结构。总共有100个随机结构，随机结构的总数是由8_ReStart.sh中ENS_NUM控制的

u_population1_1.dat, ..., u_population1_100.dat
u_population1_1.dat 代表第1代随机结构的第1个随机结构相应的位移模式
```


注意修改1_CreatEns.py中的以下参数:
```python
# nqirr = 2 代表计算2个不可约动力学矩阵
dyn = CC.Phonons.Phonons("DYNPOP", nqirr = 2) 
```

### <span style="color:yellow"> 1.3 生成qe自洽输入文件2_CreatQEInput.py
2_CreatQEInput.py是生成qe自洽输入文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R2_CreatQEInput.py，通过运行R2_CreatQEInput.py生成全部随机结构qe自洽输入文件。
注意修改2_CreatQEInput.py中的以下参数:
```python
...
    pseudo_dir = "/lustre/home/h240012/work/SSCHA/0_phon/pp" # 注意修改赝势路径
...
    ecutwfc = 80 # 注意修改截断能
    degauss = 0.02 # 注意修改degauss展宽
    nat = 216 # 不一样的q网格对应了不一样的超胞大小和原子数
    ntyp = 3  # 不一样的体系对应了不一样的元素组分
...
    Ce     140.116     Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf # 注意修改赝势名称
    Sc     44.955912   Sc.pbe-spn-kjpaw_psl.1.0.0.UPF  # 注意修改赝势名称
    H      1.00794     H.pbe-kjpaw_psl.1.0.0.UPF  # 注意修改赝势名称
...
2 2 4 0 0 0 # 注意修改k网格
...
```
R2_CreatQEInput.py获取data_ensemble_manual中全部的随结构scf_population1_*.dat, 然后在run_calculation中放置全部随机结构qe自洽输入文件。espresso_run_1.pwi, ..., espresso_run_100.pwi

注意事项：espresso_run_*.pwi 是qe的自洽输入文件，其中的 kpoints 需要测试，最大误差为0.1Ry，每个自洽用时最长为20min， 即：不需要保证k网格是收敛的参数，只要保证在0.1Ry能量误差的精度范围内用时尽可能少即可。kpoints 参数的修改不需要手动进各个 espresso_run_*.pwi 文件修改，只需要在 2_CreatQEInput.py 中修改，它会自动同步生成 R2_CreatQEInput.py. 

### <span style="color:yellow"> 1.4 生成提交任务的脚本3_Creat_Sub.py
这里是需要使用者修改参数的部分，针对你使用的机器和环境，进行修改。具体就是：header变量的内容和29行qe运行命令

如果需要计算总共100个qe的自洽，您想用20个任务，每个任务5个自洽完成，那么您直接更改3_Creat_Sub.py中的这两个值就可以。
```python
n_sub = 20
n_pw=5
```
他会帮您自动生成提交任务的脚本

### <span style="color:yellow"> 1.5 提交运算 全部随机结构的自洽.4_SubAllJobs.sh 

这里是需要使用者修改参数的部分, 具体来说：修改相应作业脚本的提交命令

## <span style="color:red"> 1.2到1.5的操作都可以用过脚本8_ReStart.sh完成
```shell
START_POP=0 #计算第0代的动力学矩阵
ENS_NUM=100 #产生100个随机结构并自洽
```

### <span style="color:yellow"> 1.6 检查自洽计算情况5_CheckFinish.py

在等待计算的过程中可以用5_CheckFinish.py检查计算情况，有多少算完的，有多少还没算的。如果全部交上去的计算任务算完，但是有部分报错或出问题，可以用5B_SubUnfinish.py重新提交。

### <span style="color:yellow"> 1.7 第j代的自洽计算完成, 收集run_calculation中第j代的所有随机结构自洽计算出的能量和受力, 以sscha可以读取的格式存储在data_ensemble_manual中 6_ParseOutput.py

前面提到了在1_CreatEns.py运行后，data_ensemble_manual已经产生了两类文件，分别是：scf_populationX_Y.dat和u_populationX_Y.dat。

在运行6_ParseOutput.py之后，data_ensemble_manual中还会产生3类文件，分别是：
```
pressures_populationX_Y.dat

forces_populationX_Y.dat

energies_supercell_populationX.dat
```

### <span style="color:yellow">  1.8 结构弛豫, 获得第j+1代的动力学矩阵 dyn_popj+1_1, dyn_popj+1_2, dyn_popj+1_3, ...... 0_Relax.py 
注意修改0_Relax.py的以下参数
```python
...
dyn = CC.Phonons.Phonons("dyn_pop11_", nqirr = 4) # "dyn_pop11_"是第j-1代的动力学矩阵 nqirr=4 动力学矩阵的个数
...
ensemble.load("data_ensemble_manual", population = 12, N = 1000) # population = 12是第j代的编号, 也就是通过弛豫会获得动力学矩阵dyn_pop12_*, N = 1000 是第j代总共有1000个结构的能量和受力要读取 (例如这里的j=12)

relax.vc_relax(target_press =150, static_bulk_modulus = 300, ensemble_loc = "data_ensemble_manual",restart_from_ens = True,start_pop = ENDPOP) # 注意修改target_press为你指定优化的压强
```

这里注意区分relax.vc_relax与relax.relax
```python
relax.relax(get_stress = True) # save the stress to compute the pressure

relax.vc_relax(target_press=PRESS, fix_volume = True) # Fixing the volume improves the convergence of the variable cell algorithm. If fix_volume=true (default False) the volume is fixed, therefore only the cell shape is relaxed.

relax.vc_relax(target_press=PRESS) # By default it is 0 (ambient pressure)

```


### <span style="color:yellow"> 1.9 提交结构弛豫任务7_SubRelax.sh 

## <span style="color:red">  1.6到1.9的操作都可以用过脚本9_ReEnd.sh完成
```shell
START_POP=0 #读取第0代的动力学矩阵
ENS_NUM=100 #收集100个随机结构并自洽
# 结合第0代的动力学矩阵 和 100个随机结构的受力和能量 迭代出第1代动力学矩阵
```
```shell
# 这些文件都没有用, 可以删掉
rm dyn_start_population*_*  dyn_end_population*_*
```

### <span style="color:yellow"> 1.10 检查动力学矩阵的受力和应力是否合理

#### 1.10.1 检查应力应变是否合理
在0_Relax.py中, 我们设置了将结构最终优化到指定压强，检查OUT*.dat中的STRESS，可以查看结构的晶格的受到的应力是否在指定压强附近, stress是否稳定在指定压强附近, 目前来看非常重要, 直接决定着最终计算出来的声子谱是否稳定. 
```shell
alias gstress='grepstress() { grep -A3 stress OUT"$1".dat | tail -n 4;  grep -A3 STRESS OUT"$1".dat | tail -n 4; }; grepstress'

# 使用方法
gstress 8
```

#### 1.10.2 检查原子受力是否趋于收敛
检查OUT*.dat中的FC，可以查看结构中原子的受力是否趋于收敛
```shell
grep FC OUT*.dat
# 如果这一代收敛的不够好，我们就进行下一代, 一直到FC gradient modulus < 0.0001

alias gfc='grepfc() { grep "FC gradient modulus" OUT"$1".dat | tail -n 1; }; grepfc'
# 使用方法: 
gstress 8
```



### <span style="color:yellow"> 1.11 至此完成了第j代动力学矩阵的迭代!!!下面介绍一些特殊情况的处理。需要注意，如果你忘记现在迭代到第几代了，下面教你一些简单的方法判断你进行到了第几代：
1. 如果8_ReStart.sh中的STRAT_POP=n1,  9_ReEnd.sh中的START_POP=n2, n1=n2且dyn_pop0_n3中的n3=n1+1=n2+1, 说明你这时已经迭代出了第n3代的动力学矩阵，下一步应该修改8_ReStart.sh中的STRAT_POP=n3然后产生随机结构做自洽了。
2. 如果8_ReStart.sh中的STRAT_POP=n1,  9_ReEnd.sh中的START_POP=n2, n1=n2+1, 说明你这时已经完成了第n1代的动力学矩阵的随机结构的自洽，下一步需要修改9_ReEnd.sh中的START_POP=n1，迭代出n1+1代的动力学矩阵了。注意代执行9_ReEnd.sh之前，请确保所有的自洽计算都完整无误。
3. OUTn4.dat中n4的编号n4始终大于9_ReEnd.sh中的START_POP=n2的编号，因为执行9_ReEnd.sh代表从第n2代的动力学矩阵中迭代出第n2+1=n4代动力学矩阵，迭代出的输出保存在OUTn4.dat中。

#### <span style="color:lightyellow"> 1.11.1 在已经完成计算的第j代, 扩大第j+1代需要自洽的结构数
当我们发现FC不再降低时，就需要扩大每一代需要自洽的结构数，此时要改以下 4个文件中的参数
```shell
# 8_ReStart.sh 中的 ENS_NUM=200
# 3_Creat_Sub.py 中的 n_sub=40 和 n_pw=5， 保证 n_sub*n_pw = ENS_NUM
# 4_SubAllJobs.sh 中的 for i in `seq 1 40` ，seq命令的终止参数 设为 n_sub=40
# 9_ReEnd.sh 中的 ENS_NUM=200
```
一般经验是:前5代, 100个结构算, 第6,7,8代, 200个结构算, 第9,10代, 500个结构算, 第11代1000个结构算.
更一般的经验是: 当FC不再下降时,就需要扩大每一代需要自洽的结构数, 一直到FC gradient modulus < 0.0001

#### <span style="color:lightyellow">  1.11.2 在已经计算完成的 第j代 中添加更多自洽结构, 将第j代的新的随机结构的自洽能量和受力 与 第j代的旧的随机结构自洽的能量和受力结合起来, 去进行结构弛豫, 生成第j+1代的dyn动力学矩阵.
1.11.2 的这种做法比 1.11.1 的做法更加复杂, 但是好处在于计算量比较小. 仅适用于第j代已经达到了收敛, 但是算出来三阶hessian矩阵有点虚频, 你想看看增加随机结构数之后, 第j代的随机结构导出的hessian矩阵的虚频的频率有没有进一步减小. 

<span style="color:lightblue"> 1.11.2.1 save.py将第j代旧的随机结构保存出来, 保存到popj_1000目录中(例如: 这里是pop12_1000)
```python
# Load the original ensemble (first population with 1000 configurations)
dyn = CC.Phonons.Phonons("dyn_pop11_", nqirr = 4)

ens = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())

# 这里i=12, 指明要加载的代数 population=12 , 指明这一代旧随机结构数 N=1000
ens.load("data_ensemble_manual", population = 12, N=1000) 

# 这里i=12, save_bin指明将第12代的结构以***.npy的形式保存, 数据保存到pop12_1000, 指明要保存的代数 population_id = 12
ens.save_bin("pop12_1000", population_id = 12) 
```

<span style="color:lightblue"> 1.11.2.2 新建一个目录2_AddTo2K(例如: 这里新目录的含义是将这一代的随机数从1K扩展到2K), 将ens.save_bin(..., ...)保存出来的目录(例如:这里要对第12代扩展随机结构数: pop12_1000), 拷贝到2_AddTo2K中. 
```shell
ls 
# 0_phon  1_100Relax

mkdir 2_AddTo2K; ls
# 0_phon  1_100Relax  2_AddTo2K
```

<span style="color:lightblue"> 1.11.2.3 将第j-1代的动力学矩阵(例如: dyn_pop11_{1..4})拷贝过来, 在此基础上, 通过脚本8_ReStart.sh执行步骤1.1~1.5
特别注意8_ReStart.sh中的以下几个参数:
```shell
START_POP=11 # 读取第j-1代的动力学矩阵.  我们用第j-1代的动力学矩阵产生更多的随机结构, 最终通过计算, 获得第j代的动力学矩阵
ENS_NUM=1000 # 这里是结构数
```

<span style="color:lightblue"> 1.11.2.4 执行9_ReEnd.sh, 因为9_ReEnd.sh可以生成R6_ParseOutput.py, 通过R6_ParseOutput.py将第j代新生成的qe格式的随机结构的能量和受力, 转化为sscha识别的能量和受力, 将结果保存在data_ensemble_manual中. 
注意这里要修改9_ReEnd.sh中的以下参数


<span style="color:lightblue"> 1.11.2.5  merg.py 合并第j代新计算出的随机结构自洽的能量和受力(存储在data_ensemble_manual中) 和 旧的随机结构的能量和受力(存储在pop12_1000), 合并好的总的随机结构的能量和受力存储在pop12_2000中.
```python
# Load the original ensemble (first population with 1000 configurations)
dyn = CC.Phonons.Phonons("dyn_pop11_", nqirr = 4) 

ens = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())

# 载入保存在pop12_1000目录下的第j代的随机结构的能量和受力(例如: 这里j=12), population_id = 12 是因为在pop12_1000目录中动力学矩阵被save.py中的ens.save_bin("pop12_1000", population_id = 12) 保存为 dyn_gen_pop12_1  dyn_gen_pop12_2  dyn_gen_pop12_3  dyn_gen_pop12_4
ens.load_bin("pop12_1000", population_id=12) 
# ens.load_bin("pop12_1000", population_id=1)  
# 如果你看到population_id=1, 是因为在pop12_1000目录中动力学矩阵被save.py中ens.save_bin("pop12_1000", population_id = 1) 的保存为 dyn_gen_pop1_1  dyn_gen_pop1_2  dyn_gen_pop1_3  dyn_gen_pop1_4

new_ensemble = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())

# 加载第j代新算出来的随机结构的能量和受力, 他们都保存在popj_1000中, 这里加载的方式不用load_bin, 是因为能量和受力被保存为*.dat的格式, 不是*.npy的格式. (例如, 这里的j=12), population = j指定加载第j代新增加的随机结构, N=1000指定结构数. 
new_ensemble.load("data_ensemble_manual", population = 12, N=1000)


# Merge the two ensembles
ens.merge(new_ensemble)

# Now ens contains the two ensembles. You can save it or directly use it for a SSCHA ncalculation
ens.save_bin("pop12_2000", population_id = 12)
```


<span style="color:lightblue"> 1.11.2.6   用 Relax_Add.py 和 R7_SubRelax.sh 计算结构弛豫, 即: 利用合并好的随机结构的能量和受力, 进行结构弛豫获得第j+1代的动力学矩阵.
注意修改Relax.py文件中以下参数:
```python
ensemble.load_bin("pop12_2000", population_id = 12) # 指明要加载的随机结构存放的目录.
relax = sscha.Relax.SSCHA(
    minim, ase_calculator = None,
    N_configs = 1000,  # 修改为2000, 2000是合并后的总的结构数
    max_pop = 1,
    save_ensemble = True,
    cluster = None,
    )
IO_freq.SetupSaving("12_freqs") 
```

```shell
# 将sbatch R7_SubRelax.sh注释了, 因为弛豫计算脚本Relax.py要修改后, 手动提交该计算任务
sbatch R7_SubRelax.sh
```

<span style="color:lightblue"> 1.11.2.6  在完成迭代计算之后可以用以下脚本来绘制收敛情况的总结图：
```shell
sscha-plot-data.py 1_freqs 2_freqs 3_freqs 4_freqs 5_freqs 6_freqs 7_freqs 8_freqs 9_freqs 10_freqs 11_freqs 12_freqs 13_freqs 14_freqs 15_freqs 16_freqs 17_freqs 18_freqs 19_freqs 20_freqs 21_freqs 22_freqs 23_freqs 24_freqs
```

## <span style="color:red"> 2. 使用受力收敛的动力学矩阵去计算三阶HESSIAN矩阵

### <span style="color:yellow"> 2.1. 当原子受力达到目标精度后(0.0001), 开始计算三阶HESSIAN矩阵

#### <span style="color:lightyellow"> 2.1.1 修改计算脚本 1_CalHess.py, 有以下 10 个参数需要修改
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

#### <span style="color:lightyellow"> 2.1.2 将1_CalHess.py提交到节点运行, 可能需要你根据自己使用的集群修改提交作业的脚本
```shell
sbatch SUB_HESS.sh
```

#### <span style="color:lightyellow"> 2.1.3 把任务提交上去之后可能爆出这样一个问题:  
```shell
ImportError: /home/h240012/soft/miniconda3/envs/sscha/lib/python3.10/site-packages/mpi4py/MPI.cpython-310-x86_64-linux-gnu.so: undefined symbol: MPI_Aint_add
```
我查了一下网上的教程, 可能是因为mpi4py相应的openmp没装好导致的, 用下面的命令就能装好, 注意千万不要用pip安装
```shell
conda install mpi4py
```

#### <span style="color:lightyellow"> 2.1.4 把任务提交上去之后还可能爆出这样一个问题:  
```shell
...
sscha/lib/pvthon3,10/site-packages/spglib/spglib.py
...
TypeError: float() argument must be a string or a real number, not 'Atom'
```
这是因为spglib版本太高与sscha不匹配, 需要降低spglib的版本
```shell
pip install spglib==1.16.5
```

#### <span style="color:lightyellow"> 2.1.5 任务计算完之后, OUT.dat中最后一行会显示DONE, 此时表明三阶HESSIAN计算完成.
可以通过以下命令查看三阶HESSIAN矩阵是否有虚频. 
```shell
for i in 1 2 3 4; do echo $i; grep freq V3_Hessian.dyn$i | head -n 10;  done
```

如果虚频太大的话, 就需要继续增加第j代的随机结构数, 尝试通过增加随机结构数, 从而提升能量和受力的精度, 降低虚频. 可能出现几种情况:

(1). 虚频越来越小, 并且收敛的速度, 趋势都可观, 说明该方法有效
(2). 虚频越来越大, 说明该结构本身就不稳定
(3). 虚频越来越小, 但是收敛的很慢, 可以尝试减小q网格, 减小不可约动力学矩阵的数量, 从而减小自洽计算的超胞, 提升收敛性.

## <span style="color:red"> 3. 使用V3_Hessian.dyn*中的结构, 计算两组动力学矩阵, 一组稀疏的, 一组稠密的, 然后用这两组动力学矩阵和V3_Hessian.dyn*插值获得最终的非谐的动力学矩阵

### <span style="color:yellow"> 3.1 重新计算稀疏的、稠密的动力学矩阵

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
python qedyn2struct.py -i V3_Hessian.dyn1 -o qe

# 然后你需要格外注意并修改scffit.in 和 scf.in 中以下部分的参数
# You have to confirm the 4 iterms, namely 
# `pp path`, 
# `ATOMIC_SPECIES`, 
# `K_POINTS in scffit.in`, 
# `K_POINTS in scf.in`.
```

如果报错，`D_S (l=2) for this symmetry operation is not orthogonal`，
这说明V3_Hessian.dyn可能弛豫出来的对称性微微扭曲，建议在那一代补充增加结构数，获得新的V3_Hessian.dyn矩阵。
```shell
task #         0
from d_matrix : error #         3
D_S (l=2) for this symmetry operation is not orthogonal
```

### <span style="color:yellow"> 3.2 回收`1.sparse`目录中稀疏q网格的动力学矩阵，`2.fine`目录中稠密q网格的动力学矩阵，`1.sscha-relax`目录中V3_Hessian动力学矩阵。
用到了`3.Inter`里面的`get_dense.sh`，`get_sparse.sh`，`get_v3_hessian.sh`三个脚本。有相关路径需要自己修改，很简单。

### <span style="color:yellow"> 3.3 在`2.interpolation/3.Inter`目录中获得插值动力学矩阵：`inter_dyn_*`

修改`2_Inter.py`的内容，并通过`SUB_INTER.sh`提交任务。
```python
···
...
# 详细说明需要修改的内容：
inter_dyn = CC.Phonons.Phonons()
dyn_sscha = CC.Phonons.Phonons("V3_Hessian.dyn",4) # 指定好V3_Hessian的名称和不可以q点数量。
dyn_coarse= CC.Phonons.Phonons("222.dyn",4) # 指定好稀疏动力学矩阵的名称和不可以q点数量。
dyn_fine  = CC.Phonons.Phonons("666.dyn",20)# 指定好稠密动力学矩阵的名称和不可以q点数量。
# 以上名称有赖于你自己的定义。前后统一即可。
inter_dyn = dyn_sscha.Interpolate(coarse_grid=[2,2,2], fine_grid=[6,6,6], support_dyn_coarse=dyn_coarse, support_dyn_fine=dyn_fine, symmetrize=True)
#dyn_222.SwapQPoints(dyn_sscha)
#dyn_222.save_qe("i")
inter_dyn.save_qe("inter_dyn_") # 插好值的动力学矩阵被命名为inter_dyn_1, inter_dyn_2, ......
···
···
```

### <span style="color:yellow"> 3.4 利用插值好的动力学矩阵计算超导

1. 首先，将`inter_dyn_*`, 降低拷贝为指定体系名称的动力学矩阵。只需要看看你之前的动力学矩阵的命名规则就行。
2. 其次在`3.Tc`目录中准备如下四个文件：`scffit.in`(自洽输入文件la2F=.true.), `scf.in`(自洽输入文件la2F=.false.), `s5_PhAssignQ.sh`(提交任务的脚本), `split_ph.in`(分q点计算的脚本)。特别注意：
   ```shell
   Electron-phonon coefficients for Nb4H14
    &inputph
    tr2_ph=1.0d-14,
    prefix='Nb4H14',
    fildvscf='Nb4H14.dv',
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


3. 从`2.fine`目录中拷贝记录了不可约q点坐标的`*.dyn0`文件到`3.Tc`目录中。
   
4. 准备分q点目录，并且在各个分q点目录中tmp文件中的`dv文件`和`patterns文件`。
   
    **这里非常重要。因为通过`ldisp=.false.`和`XQ1 XQ2 XQ3`方式指定q点，所以每一个分q点的tmp/_ph0/*.phsave/目录中patterns.*.xml都要改为tmp/_ph0/*.phsave/patterns.1.xml**

    **<span style="color:lightblue"> 以上步骤都可以通过`pre.sh`准备好, 准备好之后通过`sub.sh`提交任务。**

5. 后面的步骤就和`QE`计算超导的步骤一样了。




