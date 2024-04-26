## 流程大纲

### dyn文件改名: 222q网格计算出来的4个dyn文件分别改名为
mv Ce1Sc2H24.dyn1 dyn_pop0_1
mv Ce1Sc2H24.dyn2 dyn_pop0_2
mv Ce1Sc2H24.dyn3 dyn_pop0_3
mv Ce1Sc2H24.dyn4 dyn_pop0_4
dyn_pop0_1 代表第0代的第1个动力学矩阵，第0代也就是初始代

### 1_CreatEns.py生成随机结构
1_CreatEns.py是生成自洽文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R1_CreatEns.py, 通过运行R1_CreatEns.py生成随机结构。

R1_CreatEns.py将生成的随机结构存放在data_ensemble_manual中，分别是scf_population1_1.dat, ..., scf_population1_100.dat。
scf_population1_1.dat 代表第1代随机结构的第1个随机结构。总共有100个随机结构，随机结构的总数是由8_ReStart.sh中ENS_NUM控制的。

### 2_CreatQEInput.py生成qe自洽输入文件
2_CreatQEInput.py是生成qe自洽输入文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R2_CreatQEInput.py，通过运行R2_CreatQEInput.py生成全部随机结构qe自洽输入文件。

R2_CreatQEInput.py获取data_ensemble_manual中全部的随结构scf_population1_*.dat, 然后在run_calculation中放置全部随机结构qe自洽输入文件。espresso_run_1.pwi, ..., espresso_run_100.pwi

注意事项：espresso_run_*.pwi 是qe的自洽输入文件，其中的 kpoints 需要测试，最大误差为0.1Ry，每个自洽用时最长为20min， 即：不需要保证k网格是收敛的参数，只要保证在0.1Ry能量误差的精度范围内用时尽可能少即可。kpoints 参数的修改不需要手动进各个 espresso_run_*.pwi 文件修改，只需要在 2_CreatQEInput.py 中修改，它会自动同步生成 R2_CreatQEInput.py. 

### 3_Creat_Sub.py生成提交任务的脚本
这里是需要使用者修改参数的部分，针对你使用的机器和环境，进行修改。具体就是：header变量的内容和29行qe运行命令

如果需要计算总共100个qe的自洽，您想用20个任务，每个任务5个自洽完成，那么您直接更改3_Creat_Sub.py中的这两个值就可以。
```python
n_sub = 20
n_pw=5
```
他会帮您自动生成提交任务的脚本

### 4_SubAllJobs.sh 提交运算 全部随机结构的自洽.

这里是需要使用者修改参数的部分, 具体来说：修改相应作业脚本的提交命令

### 5_CheckFinish.py检查计算情况

在等待计算的过程中可以用5_CheckFinish.py检查计算情况，有多少算完的，有多少还没算的。如果全部交上去的计算任务算完，但是有部分报错或出问题，可以用5B_SubUnfinish.py重新提交。

### 6. 检查计算的结果是否合理

#### 6.1 检查应力应变是否合理
在0_Relax.py中, 我们设置了将结构最终优化到指定压强，检查OUT*.dat中的STRESS，可以查看结构的晶格的受到的应力是否在指定压强附近
```shell
grep -a3  "STRESS TENSOR" OUT*.dat
```

#### 6.2 检查原子受力是否趋于收敛
检查OUT*.dat中的FC，可以查看结构中原子的受力是否趋于收敛
```shell
grep FC OUT*.dat
# 如果这一代收敛的不够好，我们就进行下一代, 一直到FC gradient modulus < 0.0001
```

#### 6.3 扩大每一代需要自洽的结构数
当我们发现FC不再降低时，就需要扩大每一代需要自洽的结构数，此时要改以下 4个文件中的参数
```shell
# 8_ReStart.sh 中的 ENS_NUM=200
# 3_Creat_Sub.py 中的 n_sub=40 和 n_pw=5， 保证 n_sub*n_pw = ENS_NUM
# 4_SubAllJobs.sh 中的 for i in `seq 1 40` ，seq命令的终止参数 设为 n_sub=40
# 9_ReEnd.sh 中的 ENS_NUM=200
```
一般经验是:前5代, 100个结构算, 第6,7,8代, 200个结构算, 第9,10代, 500个结构算, 第11代1000个结构算.
更一般的经验是: 当FC不再下降时,就需要扩大每一代需要自洽的结构数, 一直到FC gradient modulus < 0.0001

### 7. 当原子受力达到目标精度后(0.0001), 开始计算三阶HESSIAN矩阵

#### 7.1 修改计算脚本1_CalHess.py, 有以下x个参数需要修改
```python
N_RANDOM = 1000 # 动力学矩阵达到收敛的那一代使用的总随机结构的数量
DYN_PREFIX =  'dyn_pop11_' # 达到收敛的那一代动力学矩阵的前一代动力学矩阵
FINAL_DYN =   'dyn_pop12_' # 达到收敛的那一代动力学矩阵
SAVE_PREFIX = 'V3_Hessian.dyn' # 计算出的Hessian矩阵的结果
NQIRR = 4 # 总的动力学矩阵的数量
Tg = 0 # 温度
T =  0 # 温度
POPULATION = 12 # 达到收敛的那一代动力学矩阵的编号, 方便从DATA_DIR="data_ensemble_manual"中读取那一代动力学矩阵计算出的受力和能量
```

#### 7.2 将1_CalHess.py提交到节点运行, 可能需要你根据自己使用的集群修改提交作业的脚本
```shell
sbatch SUB_HESS.sh
```

#### 7.3 把任务提交上去之后可能爆出这样一个问题:  
```shell
ImportError: /home/h240012/soft/miniconda3/envs/sscha/lib/python3.10/site-packages/mpi4py/MPI.cpython-310-x86_64-linux-gnu.so: undefined symbol: MPI_Aint_add
```
我查了一下网上的教程, 可能是因为mpi4py相应的openmp没装好导致的, 用下面的命令就能装好, 注意千万不要用pip安装
```shell
conda install mpi4py
```

#### 7.4 把任务提交上去之后还可能爆出这样一个问题:  
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

#### 7.5 任务计算完之后, OUT.dat中最后一行会显示DONE, 此时表明三阶HESSIAN计算完成.

### 8.  
