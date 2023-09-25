# Deep+结构预测

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

### <span style="color:yellow">  训练模型命令：
```shell
dp train input.json
```

###  <span style="color:yellow"> 远程提交作业
不在本地训练，而是在集群训练，需要准备提交作业的脚本。参考center-slurm目录下的slurm作业脚本

###  <span style="color:yellow"> 检查模型收敛性
训练完了，需要检查哪些文件去判断模型是否可以投入使用？

检查lcurve.out文件的第4、5行, 第6、7行, 分别是能量的训练和测试误差，力的训练和测试误差
```shell
#                                      能量误差     能量误差      力的误差     力的误差
#  step      rmse_val    rmse_trn    rmse_e_val  rmse_e_trn    rmse_f_val  rmse_f_trn         lr
      0      2.87e+01    5.23e+01      6.84e-01    1.55e+00      9.07e-01    1.65e+00    1.0e-03
   2000      5.59e+00    8.82e+00      1.31e-01    1.59e-01      1.77e-01    2.79e-01    1.0e-03
   4000      1.35e+01    6.88e+00      8.28e-02    2.23e-02      4.26e-01    2.18e-01    1.0e-03
   6000      4.90e+00    9.40e+00      1.46e-02    7.86e-02      1.60e-01    3.06e-01    9.4e-04
   8000      7.29e+00    6.48e+00      1.81e-01    1.31e-01      2.37e-01    2.11e-01    9.4e-04
...
...
 988000      1.09e-01    1.44e-01      1.81e-02    3.05e-02      5.50e-02    3.30e-02    1.2e-08
 990000      4.42e-02    1.93e-01      1.56e-04    3.52e-02      4.39e-02    3.54e-02    1.1e-08
 992000      6.04e-02    1.09e-01      1.09e-02    1.82e-02      3.22e-02    5.50e-02    1.1e-08
 994000      3.08e-01    1.67e-02      7.39e-02    2.73e-03      4.46e-02    1.23e-02    1.1e-08
 996000      1.08e-01    6.86e-02      1.79e-02    8.16e-03      2.36e-02    4.84e-02    1.1e-08
 998000      1.43e-01    3.28e-01      2.93e-02    7.29e-02      4.90e-02    7.90e-02    1.1e-08
1000000      5.58e-02    4.34e-02      6.83e-03    1.53e-04      3.83e-02    4.32e-02    1.0e-08

```

##  <span style="color:red"> 3. 冻结模型以及测试模型

###  <span style="color:yellow"> 冻结模型的命令

```shell
deep freeze -o graph.pb 
```

###  <span style="color:yellow"> 测试模型

```shell
# -m 指定模型
# -s 指定测试集路径
# -d 指定测试结果存储路径
dp test -m graph.pb -s ../../data/test -d result
```

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
DEEPMD INFO    # ---------------output of dp test--------------- 
DEEPMD INFO    # testing system : ../../../dp-rawdata/Lu3Li2H34
DEEPMD INFO    # number of test data : 1 
DEEPMD INFO    Energy MAE         : 3.399968e-02 eV
DEEPMD INFO    Energy RMSE        : 3.399968e-02 eV
DEEPMD INFO    Energy MAE/Natoms  : 8.717866e-04 eV
DEEPMD INFO    Energy RMSE/Natoms : 8.717866e-04 eV
DEEPMD INFO    Force  MAE         : 1.439148e-02 eV/A
DEEPMD INFO    Force  RMSE        : 2.162426e-02 eV/A
DEEPMD INFO    Virial MAE         : 3.520078e-01 eV
DEEPMD INFO    Virial RMSE        : 6.096921e-01 eV
DEEPMD INFO    Virial MAE/Natoms  : 9.025841e-03 eV
DEEPMD INFO    Virial RMSE/Natoms : 1.563313e-02 eV
DEEPMD INFO    # ----------------------------------------------- 
DEEPMD INFO    # ----------weighted average of errors----------- 
DEEPMD INFO    # number of systems : 251
DEEPMD INFO    Energy MAE         : 3.276674e-01 eV
DEEPMD INFO    Energy RMSE        : 5.042938e-01 eV
DEEPMD INFO    Energy MAE/Natoms  : 1.284612e-02 eV
DEEPMD INFO    Energy RMSE/Natoms : 2.085387e-02 eV
DEEPMD INFO    Force  MAE         : 2.561809e-02 eV/A
DEEPMD INFO    Force  RMSE        : 4.842363e-02 eV/A
DEEPMD INFO    Virial MAE         : 7.235480e-01 eV
DEEPMD INFO    Virial RMSE        : 1.919634e+00 eV
DEEPMD INFO    Virial MAE/Natoms  : 2.767808e-02 eV
DEEPMD INFO    Virial RMSE/Natoms : 8.093215e-02 eV
DEEPMD INFO    # ----------------------------------------------- 
```