# deepkit+结构预测


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
Package                 Version               Editable project location
----------------------- --------------------- -------------------------------
absl-py                 1.0.0
aiohttp                 3.8.3
aiosignal               1.2.0
ase                     3.23.0
astunparse              1.6.3
async-timeout           4.0.2
attrs                   22.1.0
backrefs                5.0.1
bcrypt                  4.2.0
blinker                 1.4
boltons                 23.0.0
bracex                  2.1.1
brotlipy                0.7.0
cachetools              4.2.2
certifi                 2023.5.7
cffi                    1.15.1
charset-normalizer      2.0.4
click                   8.0.4
cloudpickle             2.2.1
conda                   23.3.1
conda-package-handling  2.0.2
conda_package_streaming 0.7.0
contourpy               1.2.1
cryptography            39.0.1
custodian               2024.6.24
cycler                  0.12.1
dargs                   0.4.8
deepmd-kit              2.2.2
dpdata                  0.2.15
dpdispatcher            0.6.5
dpgen                   0.10.1.dev82+gdc57c9a /lustre/home/h240012/soft/dpgen
emmet-core              0.60.2
flatbuffers             1.12
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
joblib                  1.4.2
jsonpatch               1.32
jsonpointer             2.1
keras                   2.9.0
Keras-Preprocessing     1.1.2
kiwisolver              1.4.5
latexcodec              3.0.0
Markdown                3.4.1
MarkupSafe              2.1.1
matplotlib              3.9.1
mkl-fft                 1.3.1
mkl-random              1.2.2
mkl-service             2.4.0
monty                   2024.7.12
mp-api                  0.33.3
mpi4py                  3.1.3
mpmath                  1.3.0
msgpack                 1.0.8
multidict               6.0.2
networkx                3.3
numkit                  1.3.0
numpy                   1.25.0
oauthlib                3.2.0
opt-einsum              3.3.0
packaging               23.0
palettable              3.3.3
pandas                  2.2.2
paramiko                3.4.0
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
pydantic                1.10.17
PyJWT                   2.4.0
pymatgen                2022.5.18
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
six                     1.16.0
spglib                  2.5.0
sympy                   1.13.1
tabulate                0.9.0
tenacity                8.5.0
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
Werkzeug                2.2.3
wheel                   0.38.4
wrapt                   1.14.1
yarl                    1.8.1
zstandard               0.19.0

```

###  <span style="color:yellow"> deepgen的命令
```shell
nohup dpgen run param.json machine.json 1>log 2>err &
# param.json 参数设置文件
# machine.json 机器配置文件
ps -ef | grep dpgen
# -ef 显示进程具体信息
```


#我修改的部分

```python
/public/home/mayuan/miniconda3/envs/deepmd/lib/python3.10/site-packages/dpdispatcher/machines/pbs.py

#ret, stdin, stdout, stderr = self.context.block_call("qstat -x " + job_id)
ret, stdin, stdout, stderr = self.context.block_call("qstat -l " + job_id)
```


##  <span style="color:yellow"> deepgen+calypso的使用注意事项

###  <span style="color:green"> 关于record.dpgen的解释
```shell
# 第0代
00 检查准备好的输入文件 
01 train
02 analysis
03 motivation and exploration

# 第1代
10 检查准备好的输入文件 
11 train
12 analysis
13 motivation and exploration

...
...

```

如果你要从02开始续算，那么就把02后面的03 10 11 12 13 ... 都删掉了. 这表示重新做motivation and exploration

###  <span style="color:green"> 注意切换dpgen的分支

下载好dpgen的代码后，用这个切换分支
```shell
git checkout origin/devel
```

###  <span style="color:green"> 注意 param.json 一对互斥的参数

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

###  <span style="color:green"> 注意 calypso_run_model_devi.py 的运行. 

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