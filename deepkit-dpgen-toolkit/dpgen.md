# dpgen+结构预测


帖子
```shell
https://bohrium-doc.dp.tech/docs/software/CALYPSO/

```


##  <span style="color:red">  0. 安装软件，准备目录以及相应的输入文件和训练集

###  <span style="color:yellow"> 0.1 安装deepmd-kit
```shell
https://github.com/deepmodeling/deepmd-kit/releases/tag/v2.2.2
# 离线安装cpu版本
deepmd-kit-2.2.2-cpu-Linux-x86_64.sh
# 下载好之后直接
bash deepmd-kit-2.2.2-cpu-Linux-x86_64.sh
```

###  <span style="color:yellow"> 0.2 安装dpgen
```shell
https://github.com/wangzyphysics/dpgen/tree/devel

git checkout -b remotes/origin/devel
# 适用于CALYPSO+DPGEN的版本
# 下载时切记要用git clone 下载，千万不要下载安装包
pip install -e .
```

###  <span style="color:yellow"> 我这里给出一个完美的环境本版

```shell
# 检查numpy是否安装完成
python -c "import numpy, sys; sys.exit(numpy.test() is False)"
```

```shell
Package                 Version
----------------------- ---------------------
absl-py                 1.0.0
aiohttp                 3.8.3
aiosignal               1.2.0
ase                     3.23.0
astunparse              1.6.3
async-timeout           4.0.2
attrs                   22.1.0
backrefs                5.0.1
bcrypt                  4.2.0
blinker                 1.8.2
boltons                 23.0.0
bracex                  2.1.1
brotlipy                0.7.0
cachetools              4.2.2
certifi                 2023.5.7
cffi                    1.15.1
charset-normalizer      2.0.4
click                   8.1.7
cloudpickle             2.2.1
conda                   23.3.1
conda-package-handling  2.0.2
conda_package_streaming 0.7.0
contourpy               1.2.1
cryptography            39.0.1
custodian               2019.2.10
cycler                  0.12.1
dargs                   0.4.8
deepmd-kit              2.2.2
dpdata                  0.2.15
dpdispatcher            0.6.6
dpgen                   0.10.1.dev82+gdc57c9a  # 它的主程序代码位置是： deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py，跑出bug了就来看看这里。
Flask                   3.0.3
flatbuffers             2.0
fonttools               4.53.1
frozenlist              1.3.3
gast                    0.4.0
google-auth             2.6.0
google-auth-oauthlib    0.4.4
google-pasta            0.2.0
GromacsWrapper          0.9.0
grpcio                  1.48.2
h5py                    3.7.0
horovod                 0.27.0
idna                    3.4
itsdangerous            2.2.0
Jinja2                  3.1.4
joblib                  1.4.2
jsonpatch               1.32
jsonpointer             2.1
keras                   2.9.0
Keras-Preprocessing     1.1.2
kiwisolver              1.4.5
latexcodec              3.0.0
Markdown                3.4.1
MarkupSafe              2.1.1
matplotlib              3.9.2
mkl-fft                 1.3.1
mkl-random              1.2.2
mkl-service             2.4.0
monty                   2024.7.30
mpi4py                  3.1.3
mpmath                  1.3.0
multidict               6.0.2
networkx                3.3
numkit                  1.3.0
numpy                   1.25.0     #####   ----- very important -------------
oauthlib                3.2.0
opt-einsum              3.3.0
packaging               23.0
palettable              3.3.3
pandas                  2.2.2
paramiko                3.4.1
phonopy                 2.27.0     #####    ----- very important -------------
pillow                  10.4.0
pip                     23.0.1
plotly                  5.23.0
pluggy                  1.0.0
protobuf                3.20.3
psutil                  5.9.0
pyasn1                  0.4.8
pyasn1-modules          0.2.8
pybtex                  0.24.0
pycosat                 0.6.4
pycparser               2.21
PyDispatcher            2.0.7
PyJWT                   2.4.0
pymatgen                2022.5.18   #####   ----- very important -------------
PyNaCl                  1.5.0
pyOpenSSL               23.0.0
pyparsing               3.1.2
PySocks                 1.7.1
python-dateutil         2.9.0.post0
python-hostlist         1.21
pytz                    2024.1
PyYAML                  6.0
requests                2.32.3
requests-oauthlib       1.3.0
rsa                     4.7.2
ruamel.yaml             0.17.21
ruamel.yaml.clib        0.2.6
scipy                   1.14.0
setuptools              66.0.0
setuptools-scm          8.1.0
six                     1.16.0
spglib                  2.5.0
sympy                   1.13.2
tabulate                0.9.0
tenacity                9.0.0
tensorboard             2.9.0
tensorboard-data-server 0.6.1
tensorboard-plugin-wit  1.8.1
tensorflow              2.9.0
tensorflow-estimator    2.9.0
termcolor               2.1.0
tomli                   2.0.1
toolz                   0.12.0
tqdm                    4.65.0
typeguard               4.3.0
typing_extensions       4.12.2
tzdata                  2024.1
uncertainties           3.2.2
urllib3                 1.26.15
wcmatch                 8.2
Werkzeug                3.0.3
wheel                   0.38.4
wrapt                   1.14.1
yarl                    1.8.1
zstandard               0.19.0
```

###  <span style="color:yellow"> 直接复制conda环境配置dpgen

需要修改的部分：

1. `deepmd-kit/condabin/conda`的头关于python的路径修改为`realpath  deepmd-kit/bin/python`执行后的路径。
2. `deepmd-kit/bin/conda`的头关于python的路径修改为`realpath  deepmd-kit/bin/python`执行后的路径。
3. `source deepmd-kit/bin/activate`执行会报错类似`-bash: /home/h240012/soft/deepmd-kit/etc/profile.d/conda.sh: No such file or directory`。   修改`deepmd-kit/bin/activate`中关于`_CONDA_ROOT`的内容。
4. 修改`deepmd-kit/etc/profile.d/conda.sh`中关于python路径和conda路径的内容。
5. `deepmd-kit/bin/pip`的头关于python的路径修改为`realpath  deepmd-kit/bin/python`执行后的路径。
6. `deepmd-kit/bin/conda-env`的头关于python的路径修改为`realpath  deepmd-kit/bin/python`执行后的路径。这样就可以用`conda env list`查看conda环境了
7. 执行`which dp`，修改头相应的python的路径
8. 执行`which dpgen`修改头相应的python的路径
9. 执行`which dpdata`修改头相应的python的路径
10. 执行`which dpdisp`修改头相应的python的路径
11. 执行`which dp_ipi`修改头相应的python的路径
12. 执行`which dp_gmx_patch`修改头相应的python的路径
```shell
realpath  deepmd-kit/bin/python

deepmd-kit/condabin/conda
deepmd-kit/bin/conda
deepmd-kit/bin/activate
deepmd-kit/etc/profile.d/conda.sh
deepmd-kit/bin/pip
```





###  <span style="color:yellow"> dpgen+calypso运行流程

```shell
nohup dpgen run param.json machine.json 1>log 2>err &
echo $! > taskids
# param.json 参数设置文件
# machine.json 机器配置文件
ps -ef | grep dpgen
# -ef 显示进程具体信息
```

##  <span style="color:red"> 1. deepgen+calypso的使用注意事项

###  <span style="color:yellow"> 注意切换dpgen的分支

下载好dpgen的代码后，用这个切换分支
```shell
# 检查全部分支
git branch -a 

# 切换到 remotes/origin/devel
git checkout remotes/origin/devel

# 切换到 本地分支devel
git checkout devel
# Branch devel set up to track remote branch devel from origin.
# Switched to a new branch 'devel'

# 切换完了之后检查一下deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py 里面有没有关于`calypso`的关键字。确保你切换对了分支
grep calypso run.py 
```


###  <span style="color:yellow"> 关于record.dpgen的解释
```shell
# 第0代
00 make_train 检查准备好的输入文件 
01 run_train
02 post_train
03 make_model_devi
04 run_model_devi  # 即执行calypso_run_model_devi.py
05 post_model_devi # 空函数，没有任何内容
06 make_fp
07 run_fp
08 post_fp


# 第1代
10 make_train 检查准备好的输入文件 
11 run_train
12 post_train
13 make_model_devi
14 run_model_devi 
15 post_model_devi
16 make_fp
17 run_fp
18 post_fp
...
...

```

如何通过修改record.dpgen控制从指定的位置续算呢？
比如要重新运行calypso_run_model_devi.py，那么就删除03及其后面的所有，只保留00 01 02
比如要在运行完calypso_run_model_devi.py的基础上续算，那么就删除04及其后面的所有，只保留00 01 02 03

###  <span style="color:yellow">  /deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py的程序主干

```python
3791     while cont:
3792         ii += 1
3793         iter_name=make_iter_name(ii)
3794         sepline(iter_name,'=')
3795         for jj in range (numb_task) :
3796             print(jj)
3797             if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] :
3798                 continue
3799             task_name="task %02d"%jj
3800             sepline("{} {}".format(iter_name, task_name),'-')
3801             if   jj == 0 :
3802                 log_iter ("make_train", ii, jj)
3803                 make_train (ii, jdata, mdata)
3804             elif jj == 1 :
3805                 log_iter ("run_train", ii, jj)
3806                 run_train  (ii, jdata, mdata)
3807             elif jj == 2 :
3808                 log_iter ("post_train", ii, jj)
3809                 post_train (ii, jdata, mdata)
3810             elif jj == 3 :
3811                 log_iter ("make_model_devi", ii, jj)
3812                 cont = make_model_devi (ii, jdata, mdata)
3813                 if not cont :
3814                     break
3815             elif jj == 4 :
3816                 log_iter ("run_model_devi", ii, jj)
3817                 run_model_devi (ii, jdata, mdata)
3818 
3819             elif jj == 5 :
3820                 log_iter ("post_model_devi", ii, jj)
3821                 post_model_devi (ii, jdata, mdata)
3822             elif jj == 6 :
3823                 log_iter ("make_fp", ii, jj)
3824                 make_fp (ii, jdata, mdata)
3825             elif jj == 7 :
3826                 log_iter ("run_fp", ii, jj)
3827                 run_fp (ii, jdata, mdata)
3828             elif jj == 8 :
3829                 log_iter ("post_fp", ii, jj)
3830                 post_fp (ii, jdata)
3831             else :
3832                 raise RuntimeError ("unknown task %d, something wrong" % jj)
3833             record_iter (record, ii, jj)
```


###  <span style="color:yellow">  我修改的部分代码

```python
deepmd-kit/lib/python3.10/site-packages/dpdispatcher/machines/pbs.py

#ret, stdin, stdout, stderr = self.context.block_call("qstat -x " + job_id)
ret, stdin, stdout, stderr = self.context.block_call("qstat -l " + job_id)
```

```python
deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py

# standard_incar = incar_upper(Incar.from_str(incar))
standard_incar = incar_upper(Incar.from_string(incar))

# kp=Kpoints.from_str(ret)
kp=Kpoints.from_string(ret)
```


###  <span style="color:yellow"> 注意 一定要给 `2.get-dpmod/calypso_input/calypso.x`加上可执行的权限
```shell
# 不然跑不起来
chmod u+x 2.get-dpmod/calypso_input/calypso.x
```

###  <span style="color:yellow"> 注意：做变组分结构预测的时候，input.dat设置的组分变化范围太大，容易超出calypso允许产生配比的最大范围。
```shell
# 三元体系，下面的参数设置，就差不多了
# The Variation Range for each type atom
@CtrlRange
1  6
1  6
1  80
@End
```


###  <span style="color:yellow"> 注意 param.json 一对互斥的参数

```json
# 参数1
"calypso_input_path":"./calypso_input",        
"model_devi_max_iter": 15,
```

```json
# 参数2
"sys_configs": []
```

这两个参数都是可以用来控制dpgen采样的。参数1优先被读取，如果参数1没有设置，就读取参数2. 

至于参数1和参数2有什么不同，参见：https://github.com/wangzyphysics/dpgen/blob/devel/examples/run/dp-calypso-vasp/param.json

###  <span style="color:yellow"> 注意 calypso_run_model_devi.py 无法主节点运行

`iter.00000/01.model_devi/model_devi_results/`里面有个脚本`calypso_run_model_devi.py` 只能在主节点上运行, 具体是 在`iter.000000/01.model_devi/model_devi_results/` 目录下运行下面的命令:
```shell
python calypso_run_model_devi.py --all_models ../gen_stru_analy.000/graph.000.pb ../gen_stru_analy.000/graph.001.pb ../gen_stru_analy.000/graph.002.pb ../gen_stru_analy.000/graph.003.pb --type_map Ce Sc H
```

在运行过程中会产生task.000.000这个文件. 

在运行完毕后会在iter.000000/01.model_devi/model_devi_results/中产生一个Model_Devi.out的文件. 

并且在运行完毕后dpgen会自动在iter.000000/01.model_devi/record.calypso中的文件末尾会添加数字4

```shell
1 0
2 
3
4
4
```

然而在riken的主节点机器上运行这个脚本会导致一堆core-python文件的爆出. 所以最好还是把它放到远程机器上运行. 怎么操作呢？

**第1种情况：挂在后台的dpgen进程中断了**

在执行上述`python calypso_run_model_devi.py ....`命令前, 先在iter.000000/01.model_devi/record.calypso中添加数字4, 然后手动执行命令sbatch sub_devi.sh, 具体sub_devi.sh的内容已经写在各个机器相应的目录中了, 然后执行` nohup dpgen run param.json machine.json 1>log 2>err &
`

**第2种情况：挂在后台的dpgen进程正常运行**

直接手动执行命令sbatch sub_devi.sh，dpgen会自动在iter.000000/01.model_devi/record.calypso中添加数字4

最后dp官网还给出了一种命令行参数的方式获得Model_Devi.out, 详细参考：https://docs.deepmodeling.com/projects/deepmd/en/r2/test/model-deviation.html
```shell
dp model-devi -m graph.000.pb graph.001.pb graph.002.pb graph.003.pb -s ./data -o Model_Devi.out
```
<span style="color:lightblue"> **注意: 当你发现`task.*.000`中没有`model_devi.out`时，可以通过运行`dp model-devi ...`或者`python calypso_run_model_devi.py ...`获得, 然后通过执行`merge_model_devi.py`将所有的`task.*.000`中的`model_devi.out`合并为`Model_Devi.out`。具体步骤如下：**

<span style="color:lightblue"> **1. 进入指定的目录，执行`check_model_devi.sh`检查哪些目录中缺少`model_devi.out`
```shell
cd 2.getdp-mod/iter.00000*/01.model_devi/model_devi_results
./check_model_devi.sh 000 820
```
<span style="color:lightblue"> **2. 备份所有的task.*.000，这一步非常重要，因为每次运行`python calypso_run_model_devi.py ...`，其它task目录总会莫名其妙少了`model_devi.out`，为了防止你不停的丢失文件，你可以备份它们，如果真丢失了`model_devi.out`可以直接拷贝**
```shell
mkdir backup
cp task* backup
```
<span style="color:lightblue"> **3. 对需要重新运行`python calypso_run_model_devi.py ...`的结构，执行下面的命令。 另外：`python calypso_run_model_devi.py ...`运行输出的model_devi.out比`dp model-devi ...`运行输出的结果多一列`min_dis`, 相应的`python calypso_run_model_devi.py ...`运行输出的Model_Devi.out也多一列`min_dis`**

```shell
./retry_model_devi.sh
```

<span style="color:lightblue"> **4. 合并`model_devi.out`到`Model_Devi.out`. 执行下面的命令:**

```shell
python merge_model_devi.py 
```



###  <span style="color:yellow"> 注意设置压强

关于压强的参数有3个地方需要设置：
```shell
# 第1个地方
0.vasp-dataset/INCAR_{1..5}
...
PSTRESS = 2000
...

# 第2个地方
2.get-dpmod/calypso_input/input.dat
...
PSTRESS = 2000
...

# 第3个地方
2.get-dpmod/vasp_input/INCAR
...
PSTRESS = 2000
...
```

###  <span style="color:yellow"> 关于dpdata的LabeledSystem的使用小贴士：
1. `/home/h240012/soft/deepmd-kit/lib/python3.10/site-packages/dpdata/system.py:1109`是LabeledSystem的位置。
2. 索引读取后的数据的性质时，有：`energies`, `forces`, `virials`三项可以通过字典的方式索引。

###  <span style="color:yellow"> dp test 无法主节点运行，报错`Aborted (core dumped)`

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
将上述4个参数设置好后，calypso_run_model_devi.py就也可以在后台运行了。

### <span style="color:yellow"> 发现fp时候或者model_devi的时候，不报任何错，但是task任务没有跑起来。

<span style="color:lightblue">  **原因1：可能是因为节点损坏，任务提交上去后立刻掉。可以修改machine.json里面的`"custom_flags": ["SBATCH --exclude=node21"]`, 将坏节点排除出去。**

<span style="color:lightblue">  **原因2：可能是因为在写machine.json的fp部分的command参数时，加入了很多带#的内容, 直接把后面的shell命令注释了，所以最好还是不要加任何#**

### <span style="color:yellow"> 发现fp时候, 有些fp迟迟算不完。你需要进到目录里面看一下，但是你要找不到是相应的哪个task目录。
```shell
第一步：找到相应任务号
第二步：vi log
第三步：搜索相应任务号，并找到对于的哈希值x
第四步：在 2.getdp-mod/execute/fp/哈希值目录y 中找到对应的哈希值x作业脚本：x.sub.run。
里面就写着task目录。

例如：
cd $REMOTE_ROOT
cd task.100.000000
test $? -ne 0 && exit 1
if [ ! -f db5409944b5165c01038009c21b4dee6722a9365_task_tag_finished ] ;then
...
...

大功告成。
```

### <span style="color:yellow"> param.json 里面有一块参数（如下）会影响numb_steps等参数的设置。如果你要从头跑dpgen，最好把它们删掉了。

```json
{
  ...
    "training_init_model":              false,
    "training_reuse_iter":              1,
    "training_reuse_old_ratio":         0.9,
    "training_reuse_start_lr":          1e-4,
    "training_reuse_stop_batch":        1000000,
    "training_reuse_start_pref_e":      0.2,
    "training_reuse_start_pref_f":      100,
  ...
}
```

## <span style="color:red"> 2. 如何判断势训练的优劣

### <span style="color:yellow"> 检查candidate的数量

dpgen.log里面`iter.*task 06`记录了每一代相应的candidate的数量。

其中`"model_devi_f_trust_lo": 0.5`和`"model_devi_f_trust_hi": 0.9`控制着被选为candidate的数量，低于`"model_devi_f_trust_lo"`则认为是accurate，高于`model_devi_f_trust_hi`则认为是failed




### <span style="color:yellow"> 绘制The model deviation distribution
`plot_ModelDevi.py`可以用于绘制Model_devi.out中力的统计分布数。

它有两种运行模式：
```shell
You can run it by: 

python plot_corelation.py ../../../../1.dp-data/2.testset/Ce1Sc2H22/ frozen_model.pb 

or 

python plot_corelation.py ../../../../1.dp-data/2.testset frozen_model.pb
```

### <span style="color:yellow"> dptest检查数据集

```shell
dp test -m frozen_model.pb -s ../../../../1.dp-data/2.testset/  -d result > result.log 2>&1
```

主要集中看result.log的最后关于`weighted average of errors`的部分的`RMSE/Natoms`，例如：
```shell
DEEPMD INFO    # ----------weighted average of errors-----------
DEEPMD INFO    # number of systems : 29
DEEPMD INFO    Energy MAE         : 4.879041e-01 eV
DEEPMD INFO    Energy RMSE        : 7.563694e-01 eV  
DEEPMD INFO    Energy MAE/Natoms  : 8.717066e-03 eV
DEEPMD INFO    Energy RMSE/Natoms : 1.380319e-02 eV   # -------very important------
DEEPMD INFO    Force  MAE         : 1.690654e-01 eV/A
DEEPMD INFO    Force  RMSE        : 2.512424e-01 eV/A # -------very important------
DEEPMD INFO    Virial MAE         : 8.679276e+00 eV
DEEPMD INFO    Virial RMSE        : 2.049670e+01 eV
DEEPMD INFO    Virial MAE/Natoms  : 1.501872e-01 eV    
DEEPMD INFO    Virial RMSE/Natoms : 3.432773e-01 eV   # -------very important------
DEEPMD INFO    # -----------------------------------------------

```


## <span style="color:red"> 报错解读

### <span style="color:yellow"> 做fp的时候遇到如下报错，多半是因为intel的oneapi编译器在source的时候，没有加上`--foce`。在`machine.json`中的fp的部分的`command`部分写完整的source命令即可。例如：`source /public/home/liuhanyu/intel/oneapi/setvars.sh --force`.

```shell
Traceback (most recent call last):
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpdispatcher/submission.py", line 358, in handle_unexpected_submission_state
    job.handle_unexpected_job_state()
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpdispatcher/submission.py", line 862, in handle_unexpected_job_state
    raise RuntimeError(err_msg)
RuntimeError: job:734273b19c6e20d77b77e4e34ece5bb2e0a08668 530904 failed 3 times.
Possible remote error message: ==> /public/home/liuhanyu/workplace/mayuan/40.La-Sc-H/2.getdp-mod/execute/fp/9e16ce64ccdfa65bee93be36838740e338edd476/task.010.000010/fp.log <==
red on exiting setvars.sh.
  
 
:: WARNING: setvars.sh has already been run. Skipping re-execution.
   To force a re-execution of setvars.sh, use the '--force' option.
   Using '--force' can result in excessive use of your environment variables.
  
usage: source setvars.sh [--force] [--config=file] [--help] [...]
  --force        Force setvars.sh to re-run, doing so may overload environment.
  --config=file  Customize env vars using a setvars.sh configuration file.
  --help         Display this help message and exit.
  ...            Additional args are passed to individual env/vars.sh scripts
                 and should follow this script's arguments.
  
  Some POSIX shells do not accept command-line options. In that case, you can pass
  command-line options via the SETVARS_ARGS environment variable. For example:
  
  $ SETVARS_ARGS="ia32 --config=config.txt" ; export SETVARS_ARGS
  $ . path/to/setvars.sh
  
  The SETVARS_ARGS environment variable is cleared on exiting setvars.sh.
  


The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/bin/dpgen", line 8, in <module>
    sys.exit(main())
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpgen/main.py", line 185, in main
    args.func(args)
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py", line 3944, in gen_run
    run_iter (args.PARAM, args.MACHINE)
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py", line 3827, in run_iter
    run_fp (ii, jdata, mdata)
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py", line 3187, in run_fp
    run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _vasp_check_fin,
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpgen/generator/run.py", line 3166, in run_fp_inner
    submission.run_submission()
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpdispatcher/submission.py", line 231, in run_submission
    self.handle_unexpected_submission_state()
  File "/public/home/liuhanyu/workplace/mayuan/software/deepmd-kit/lib/python3.10/site-packages/dpdispatcher/submission.py", line 362, in handle_unexpected_submission_state
    raise RuntimeError(
RuntimeError: Meet errors will handle unexpected submission state.
Debug information: remote_root==/public/home/liuhanyu/workplace/mayuan/40.La-Sc-H/2.getdp-mod/execute/fp/9e16ce64ccdfa65bee93be36838740e338edd476.
Debug information: submission_hash==9e16ce64ccdfa65bee93be36838740e338edd476.
Please check error messages above and in remote_root. The submission information is saved in /public/home/liuhanyu/.dpdispatcher/submission/9e16ce64ccdfa65bee93be36838740e338edd476.json.
For furthur actions, run the following command with proper flags: dpdisp submission 9e16ce64ccdfa65bee93be36838740e338edd476

```

