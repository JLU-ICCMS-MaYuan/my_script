# 如何使用USPEX
## 1. 参考网址
1. http://www.pwmat.com:8080/upload/module/pdf/PWmatUSPEX.pdf
2. https://uspex-team.org/zh/uspex/overview
3. https://uspex-team.org/online_utilities/uspex_manual_release/ChineseVersion/uspex_manual_chinese_V10.2/index.html
4. https://uspex-team.org/online_utilities/uspex_manual_release/ChineseVersion/uspex_manual_chinese_V10.2/sect0018.html (变组分结构预测)
5. https://uspex-team.org/online_utilities/uspex_manual_release/ChineseVersion/uspex_manual_chinese_V10.2/faq_compositions.html (如何控制成分)
6. https://www.bilibili.com/read/cv8674987/
## 2. 准备好输入文件后，提交任务
```shell
nohup ./submit.sh > submit.log 2>&1 &
```

## 3. 如何准备种子文件
参考网址：https://uspex-team.org/online_utilities/uspex_manual_release/ChineseVersion/uspex_manual_chinese_V10.2/faq_seeds.html

1. 生成的结构的种子文件
2. Seeds 文件夹，这个文件夹就是 USPEX 计算时，需要种子文件存放的地方
3. 在VASP5格式中的级联POSCAR文件的格式如下，为下一代结构创建 Seeds/POSCARS文件，或为一次USPEX计算的特殊代创建Seeds/POSCARS_gen(gen指代号)。 在此文件名中不要遗忘字母“S”
4. 在计算中的任何时候都可以加入种子。每一代开始时USPEX将尝试读取种子文件的两个类型。 相应的信息将被记录在results1/Seeds_history中，种子文件（POSCARS 或 POSCARS_gen）将为被作为POSCARS_gen保存在种子Seeds/文件中。
5. 当种子加入时，我们建议用户检查results1/Seeds_history和警告文件。 如果你的种子有问题时，可能会有一个警告信息“Meet a problem when reading Seeds - ...”。 当一个错误出现在种子文件，比如是缺失了一行，错误出现后的结构将不会被添加进去。
```shell
EA33 2.69006 5.50602 4.82874 55.2408 73.8275 60.7535 no SG
1.0
2.690100  0.000000  0.000000
2.690100  4.804100  0.000000
1.344900  2.402100  3.967100
Mg Al O
1  2  4
Direct
0.799190  0.567840  0.859590
0.793520  0.230950  0.544750
0.793540  0.916090  0.174450
0.050972  0.816060  0.859610
0.172230  0.194810  0.859600
0.438250  0.655170  0.406880
0.438230  0.202440  0.312330
EA34 7.61073 2.85726 2.85725 60.0001 79.1809 79.1805 no SG
1.0
7.610700  0.000000  0.000000
0.536350  2.806500  0.000000
0.536330  1.352000  2.459300
Mg Al O
1  2  4
Direct
0.708910  0.507440  0.068339
0.374050  0.285730  0.846630
0.023663  0.069185  0.630090
0.889560  0.780560  0.341460
0.350470  0.626920  0.187820
0.597290  0.211310  0.772210
0.116440  0.371590  0.932500
......
......
......
```

## 4. Submission目录中的注意事项
1. submitJob_local.py 文件中有部分代码需要根据集群提交任务的系统不同，需要相应的修改

```shell
# slurm 系统
    RUN_FILENAME = 'myrun'
    JOB_NAME = 'optrun-{}'.format(index)
 
    # Step 1
    myrun_content = '''#!/bin/sh
#!/bin/bash
#SBATCH  --job-name={}
#SBATCH  --output=opt.out
#SBATCH  --error=opt.err                                  
#SBATCH  --partition=lhy
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1   
 
source /work/env/intel2018
  
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0

{}
'''.format(JOB_NAME, commandExecutable)  
......
    output = check_output('sbatch {}'.format(RUN_FILENAME), shell=True, universal_newlines=True)                                                                               
......
    jobNumber = int(output.split(' ')[3])
......
```

```shell
# pbs 系统
    RUN_FILENAME = 'myrun'
    JOB_NAME = 'USPEX-{}'.format(index)

    # Step 1
    myrun_content = '''#!/bin/bash
#PBS -N {}
#PBS -q liuhy
#PBS -l nodes=1:ppn=28,walltime=4:00:00
#PBS -j oe
#PBS -V
source /public/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited
cd $PBS_O_WORKDIR
killall -9 vasp_std

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

{}
'''.format(JOB_NAME, commandExecutable)

......
    output = check_output('qsub {}'.format(RUN_FILENAME), shell=True, universal_newlines=True)

......
    jobNumber = int(output.split('.')[0])

......
```

## 报错处理

### Error using python_uspex (line 91) System error: Traceback (most recent call last):   File "~/USPEX/application/archive/src/FunctionFolder/random_topology.py", line 3, in <module>  from randomTopology import generate_structure ModuleNotFoundError: No module named 'randomTopology' Command executed: python3 -W ignore /public/home/liuhanyu/workplace/mayuan/software/USPEX/application/archive/src/FunctionFolder/random_topology.py 0 168.27 NONE 3 11   3   0 Error in random_topology (line 10) Error in Random_Init_301 (line 53) Error in initialize_POP_STRUC_301 (line 129) Error in Initialize (line 46) Error in Start (line 46) Error in USPEX (line 39) MATLAB:python:ExecutionError  

```shell
安装python3.8.13版本，更高版本都会报错~

conda create -n uspex python=3.8.13

conda activate uspex

```
然后按照下面的库安装. 如果你不想安装，我还有一个完整的conda环境，里面的所有库都已经配置好了名为`uspex-condaenv`，只需要解压到miniconda的env目录下即可。
```shell
Package             Version
------------------- ---------
certifi             2024.8.30
importlib_resources 6.4.5
mkl-fft             1.3.0
mkl-random          1.1.1
mkl-service         2.3.0
numpy               1.19.2
pip                 24.2
scipy               1.6.2
setuptools          75.1.0
six                 1.16.0
spglib              1.16.2
TBB                 0.2
typing_extensions   4.12.2
wheel               0.37.1
zipp                3.20.2
```