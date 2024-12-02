# deepkit+结构预测

帖子
```shell
https://bohrium-doc.dp.tech/docs/software/CALYPSO/

```

##  <span style="color:red">  0. 准备目录

```shell
mkdir -p 1_scf_datasets/1_traindatas
mkdir -p 1_scf_datasets/2_testdatas
mkdir -p 1_scf_datasets/3_createdatas

mkdir -p 2_dp_rawdata/1_traindatas
mkdir -p 2_dp_rawdata/2_testdatas

mkdir    3_module

```

## <span style="color:red">  1. 准备数据

###  <span style="color:yellow"> 准备的数据需要保证的收敛精度
```shell
total energy < 1meV/atom
force < 10meV/A
virial tensors < 1meV/atom
```

###  <span style="color:yellow"> 自洽计算
测试encut, kspacing, 精度参考上图

如果要提取能量，力，维里力的信息，可以参考extract_info.py

如果是三元体系，例如Lu-Li-H，最好准备Lu, Li, H, Lu-H, Li-H, Lu-Li的结构预测数据
四元同理

### <span style="color:yellow"> 转化数据

将计算好的vasp自洽数据通过脚本convertVasp2DpRaw.py转化为dp可以读取的数据

##  <span style="color:red"> 2. 训练模型

### <span style="color:yellow"> 输入文件input.json
开始训练模型之前，需要准备好
```shell
# 训练步数和衰减步数的关系是：
numb_steps / 200 = decay_steps
# **通常训练numb_steps的量级在百万级**
```


### <span style="color:yellow">  训练模型命令：
```shell
# deepkit的命令
# 最普通的开始训练
dp train input.json > dp.log 2>&1


# 下面的内容都是关于如何需算deepkit
# Initialize the training from the frozen model.
dp train input.json --init-frz-model graph.pb > dp.log 2>&1

# initializes the model training with an existing model that is stored in the path prefix of checkpoint files model.ckpt, 
# the network architectures should match.
dp train input.json --init-model model.ckpt > dp.log 2>&1

# 续算的话，指定续算的文件即可， Restart the training from the provided checkpoint
dp train --restart model.ckpt input.json > dp.log 2>&1
```
更多关于续算的命令参考：https://docs.deepmodeling.com/projects/deepmd/en/r2/train/training-advanced.html

##  <span style="color:red"> 3. 冻结模型以及测试模型

###  <span style="color:yellow"> 冻结模型的命令

```shell
dp freeze -o graph.pb 
```

###  <span style="color:yellow"> 测试模型

```shell
# -m 指定模型
# -s 指定测试集路径
# -d 指定测试结果存储路径
dp test -m graph.pb -s ../../data/test -d result
```

结果如下
```shell
DEEPMD INFO    # ---------------output of dp test--------------- 
DEEPMD INFO    # testing system : ../../../dp-rawdata/Lu4Li2H33
DEEPMD INFO    # number of test data : 1 
DEEPMD INFO    Energy MAE         : 2.834020e-02 eV
DEEPMD INFO    Energy RMSE        : 2.834020e-02 eV
DEEPMD INFO    Energy MAE/Natoms  : 7.266718e-04 eV
DEEPMD INFO    Energy RMSE/Natoms : 7.266718e-04 eV
DEEPMD INFO    Force  MAE         : 3.306630e-02 eV/A
DEEPMD INFO    Force  RMSE        : 5.531320e-02 eV/A
DEEPMD INFO    Virial MAE         : 6.104632e-01 eV
DEEPMD INFO    Virial RMSE        : 1.586536e+00 eV
DEEPMD INFO    Virial MAE/Natoms  : 1.565290e-02 eV
DEEPMD INFO    Virial RMSE/Natoms : 4.068042e-02 eV
DEEPMD INFO    # ----------------------------------------------- 
.....
.....
```

预测数据与原始数据之间的相关性可视化, 可以通过执行脚本plotrelation.py获得
```shell
python plotrelation.py
```

#### <span style="color:yellow"> 模型精度问题
使用MLP做分子动力学模拟时，精度的标准是什么：
1. 能量误差 < 1meV/atom
2. 力误差0.05 - 2 eV/Angstrom
3. 力误差降低到训练数据集中力的平均值的1/10以下，可以认为模型训练有精度。
4. 多时候要根据研究的问题，检验模型在实际物理性质上的预测精度。

当精度没有达到预期的时候，应该怎么办：
1. DFT的数据精度不够：garbage in, garbage out
2. 缺失数据
3. 数据中原子距离过近



###  <span style="color:yellow"> 检查模型收敛性遇到的问题
训练完了，需要检查哪些文件去判断模型是否可以投入使用？

检查lcurve.out文件的第4、5行, 第6、7行, 分别是能量的训练和测试误差，力的训练和测试误差
```shell
head -n 2 lcurve.out && tail -n 2 lcurve. out
#                                      能量误差     能量误差      力的误差     力的误差
  ****  step      rmse_val    rmse_trn    rmse_e_val  rmse_e_trn    rmse_f_val  rmse_f_trn         lr
      0      2.87e+01    5.23e+01      6.84e-01    1.55e+00      9.07e-01    1.65e+00    1.0e-03
 998000      1.43e-01    3.28e-01      2.93e-02    7.29e-02      4.90e-02    7.90e-02    1.1e-08
1000000      5.58e-02    4.34e-02      6.83e-03    1.53e-04      3.83e-02    4.32e-02    1.0e-08
```
通过执行脚本 plotlcurve.py可以可视化误差的收敛性
```shell
python plotlcurve.py
```


###  <span style="color:yellow"> dp test 无法主节点运行dp test，报错`Aborted (core dumped)`

可能是因为这些机器限制某一个任务可以使用的最大进程数，并不是不能运行，只要把下面这三个参数调小一点就可以跑起来了
更多详细的设置可以参考：`https://docs.deepmodeling.com/projects/deepmd/en/r2/troubleshooting/howtoset_num_nodes.html`



```shell
# ---- dpgen ----
# 默认DP_INFER_BATCH_SIZE=1024, 如果设置的过大，会导致`Aborted (core dumped)`
export DP_INFER_BATCH_SIZE=2048

export OMP_NUM_THREADS=3 # 设置用于 OpenMP 并行计算的线程数。
export TF_INTER_OP_PARALLELISM_THREADS=3 # 设置 TensorFlow 内部操作（Intra-Op）并行计算的线程数。
export TF_INTER_OP_PARALLELISM_THREADS=2 # 设置 TensorFlow 操作之间（Inter-Op）并行计算的线程数。
```

