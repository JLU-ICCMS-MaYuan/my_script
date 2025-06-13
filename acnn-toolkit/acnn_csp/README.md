# ACNN

## 官方文档：

https://bonjour221.github.io/notes.github.io/external/external-utilities/

##  <span style="font-size: 30px; color: lightgreen;"> 安装教程

可能需要的安装包准备： gcc(9以上版本)， ase，  pymatgen，  cmake (3以上版本)， qhull  openblas airss  xmgrace  libtorch torchdemo(ACNN) lammps(安装到torchdemo的interface里)

1. https://www.mtg.msm.cam.ac.uk/Codes/AIRSS;
2. https://plasma-gate.weizmann.ac.il/pub/grace/src/; 
3. https://www.openblas.net

安装提醒：ase, pymatgen, qhull, airss 都可以使用 conda install 安装，具体命令为：
```shell
conda install -c conda-forge ase, pymatgen, qhull, airss
```
离线安装方法为：在有线小机器上conda安装，使用 conda pack 打包，上传即可
* 如果ase 无法使用conda安装，单独使用pip install 安装即可
 
###  <span style="font-size: 30px; color: red;">  编译流程：

####  <span style="font-size: 25px; color: blue;"> 1. 编译openblas(安装完成后不需要添加export路径)

网络帖子: https://zhuanlan.zhihu.com/p/631348362

```shell
make -j8 USE_OPENMP=1
make install PREFIX=/...
```

如果在make的时候报错，比如缺少：

就写一个叫env_gcc.sh的脚本，内容包含：如何source你的gcc以及相应的lib库, 这里给出一个例子
```shell
source /public/env/gcc-9.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/gcc-9.2.0/lib64 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/gmp-6.1.2/lib 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mpfr-4.0.1/lib 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mpc-1.1.0/lib 
```

如果报错：`error while loading shared libraries: libmpfr.so.6: cannot open shared object file: No such file or directory`
要么你按照这个帖子自己安装一个gcc.9.2: https://www.zhihu.com/people/ma-yuan-94-83/posts. 要么你自己安装一个mpfr.

####  <span style="font-size: 25px; color: blue;"> 2. 安装airss(安装完成后必须添加export路径)


<span style="font-size: 20px; color: lightblue;"> 1. 上传压缩包airss-v0.9.4.tgz并解压
```shell
$ tar -xzvf airss-v0.9.4.tgz

$ ls
bin  CONTRIBUTORS  examples  external  include  lib  LICENCE  makefile  README.md  src  test  VERSION.md
```

<span style="font-size: 20px; color: lightblue;"> 2. 手动下载spglib、symbol的压缩包到external文件夹中

下载地址在external/(spglib or symbol)/makefile里有，

压缩包的名字要求spglib的压缩包名称是：v1.14.1.zip，symbol的压缩包名称是：symmol.zip

<span style="font-size: 20px; color: lightblue;"> 3. 安装之前拷贝好手动copy三个静态库（libblas.a liblapack.a libsymspg.a 在three_lib.zip里）到lib文件夹中

```shell
make

# 必须install， 不然会报出没有cable这个命令
make install
```
<span style="font-size: 20px; color: lightblue;">  意外

1. 如果没有网络, 手动下载spglib、symbol的压缩包并拷贝到external文件夹中, 注意对spglib、symbol的版本有着极其严苛的要求
你最好在一个有网络的机器上, 执行make，搞到这两个安装包，然后把它拷贝到你要安装的机器上
2. 有时候安装报错是因为没有找到关于lpack和blas的两个静态库。我也不知道怎么安装，但是聪明的师弟已经搞好了，发给我了。
你只需要把他们拷贝到lib目录下即可。（libblas.a  liblapack.a  libsymspg.a）


####  <span style="font-size: 25px; color: blue;"> 3. 安装grace(安装完成后必须添加export路径)


上传grace-latest.tar.gz到服务器，解压，进入文件夹，执行：

```shell
source ~/bin/env_gcc.9.2.0.sh
xmgrace: ./configure --prefix=... ; make -j8 ; make install
```

####  <span style="font-size: 25px; color: blue;"> 4. 安装libtorch(安装完成后不需要添加export路径)

PyTorch官网下载: https://pytorch.org/
    
![Download libtorch](picture/libtorch.jpg)
        
解压后即可。在后续安装torchdemo的时候，会用到libtorch的路径。

####  <span style="font-size: 25px; color: blue;"> 5. 安装cmake3(安装完成后必须添加export路径)

下载好安装包之后直接将其目录下的bin目录加入环境变了就可以使用

####  <span style="font-size: 25px; color: blue;"> 6. 安装acnn(安装完成后必须添加export路径)

<span style="font-size: 20px; color: lightblue;"> 1. 解压进入torchdemo

<span style="font-size: 20px; color: lightblue;"> 2. 修改prefix.cmake
```shell
# 这一部分取消注释并且修改相应的路径
set(CMAKE_CXX_COMPILER /work/software/gcc-9.2.0/bin/g++)                 # required
set(Torch_DIR /work/home/mayuan/software/libtorch)                       # required
set(OpenBLAS_DIR /work/home/mayuan/software/OpenBLAS-0.3.28/anzhuang)    # required
```

**同时特别注意，一定要将其它所有的部分都注释，不然会报错：**
```shell
-- Release
CMake Warning (dev) at CMakeLists.txt:29 (find_package):
  Policy CMP0146 is not set: The FindCUDA module is removed.  Run "cmake
  --help-policy CMP0146" for policy details.  Use the cmake_policy command to
  set the policy and suppress this warning.

This warning is for project developers.  Use -Wno-dev to suppress it.

CMake Error at CMakeLists.txt:44 (find_package):
  By not providing "FindOpenBLAS.cmake" in CMAKE_MODULE_PATH this project has
  asked CMake to find a package configuration file provided by "OpenBLAS",
  but CMake did not find one.

  Could not find a package configuration file provided by "OpenBLAS" with any
  of the following names:

    OpenBLASConfig.cmake
    openblas-config.cmake

  Add the installation prefix of "OpenBLAS" to CMAKE_PREFIX_PATH or set
  "OpenBLAS_DIR" to a directory containing one of the above files.  If
  "OpenBLAS" provides a separate development package or SDK, be sure it has
  been installed.

```

<span style="font-size: 20px; color: lightblue;"> 3. 安装

(除前五行需要修改，后面的全注释)，创建文件夹build，进入后cmake
```shell
cmake -B build
cmake --build build --target acnn
```


####  <span style="font-size: 25px; color: blue;"> 7. 安装acnn和lmp_mpi的接口(安装完成后必须添加export路径)

<span style="font-size: 20px; color: lightblue;"> 1. 进入torchdemo的interface/lammps

<span style="font-size: 20px; color: lightblue;"> 2. 安装

自己下载个压缩包(lammps-2Aug2023.tar.gz)，更改build_lammps_interface.sh中的压缩文件名之后，

激活intel进行编译（sh build_lammps_interface.sh build 核数），lmp_mpi在文件夹lammps-acnn/build里
```shell
source /work/env/intel2024
source ~/bin/env_gcc.sh 
source /work/env/cmake-3.23

sh build_lammps_interface.sh build 8
```

####  <span style="font-size: 25px; color: blue;"> 8. 安装qhull(安装完成后必须添加export路径)
直接下载好代码，make安装之后添加bin目录到bashrc

####  <span style="font-size: 25px; color: blue;"> 9. 范例：所有添加了PATH路径的代码
```shell
export PATH=$PATH:/work/home/mayuan/software/airss/bin
export PATH=$PATH:/work/home/mayuan/software/grace-5.1.25/anzhuang/grace/bin
export PATH=$PATH:/work/home/mayuan/software/torchdemo-v3/build
export PATH=$PATH:/work/home/mayuan/software/torchdemo-v3/interface/lammps/lammps-acnn/build
export PATH=$PATH:/work/home/mayuan/software/qhull-2020.2/bin

export PATH=$PATH:/work/home/mayuan/software/torchdemo-v3-old/interface/airss

export PATH=/work/home/mayuan/software/torchdemo-v3/interface/bfgs/build:$PATH
```

####  <span style="font-size: 25px; color: blue;"> 10. 使用acnn+CSP

<span style="font-size: 20px; color: lightblue;"> 1. 创建任务：
```shell
# 特别注意这个距离是用来筛选优化后不合理的结构的，设置的参数标准可以宽松一点，比如这里H-H只要保证大于0.75即可，没有两个氢原子贴的特别近，即可。
# 这个距离不同于airss里面设置的原子间距离，airss里面设置的原子间距离用于生成结构的，最好保证生成的结构更加贴近于优化后的结构，设置的参数标准可以严格一点，比如这里H-H保证生成的是原子氢，最好大于1.0。

acnn_deploy -p 200 -s CeScH -b Ce-Ce=1.88853,Ce-Sc=1.870015,Ce-H=1.351595,Sc-Sc=1.8515,Sc-H=1.33308,H-H=0.75 -n public
```

<span style="font-size: 20px; color: lightblue;"> 2. 修改参数：

```shell
# 修改ICNAR，ENCUT和KSPACING最重要
vi DFT/dyn_vasp_in  

# 修改编译环境和提交任务脚本的头文件以及vasp_std的路径
vi DFT/sub.sh

# 在DFT目录下准备POTCAR: 准备POTCAR-元素 (用cat复制，不要用cp)

# 准备airss生成文件的脚本CeScH.cell, 特别注意名字必须是与`acnn_deploy -s`指定的名字一致。
vi RSS/Base/CeScH.cell
#VARVOL=8
#MAXTIME=0.1
##NOCOMPACT

#### FORMULA ####
##SPECIES=Ce%NUM=1-2,Sc%NUM=1-2,H%NUM=3-30
#SPECIES=Ce,Sc,H
#FOCUS=3
#NATOM=3-32
#MINSEP=0.75 Ce-Ce=3.5 Ce-Sc=3.0 Ce-H=2.0 Sc-Sc=3.0 Sc-H=1.7 H-H=0.9

#### SYMMETRY ####
#NFORM=1
#SYMMOPS=1-48
##SYMMNO=2-230

# 修改机器学习势提交脚本
vi POT/sub.sh

# 修改机器学习训练圈数，即：POT/tr的nbatch
vi POT/tr 

# 修改RELAX控制提交作业的脚本，该信息存储在TASK变量中, 特别：激活的环境、编译器路径，每一代优化多少个结构
# 使用lammps优化结构
vi RELAX/dyn_batch_relax
range="1 500"

# 使用ares优化结构
vi RELAX/dyn_batch_relax_bfgs
# parallel hierarchy
frame="500" # 总共提取500个结构
group="100" # 每100个为一组，相应的n_groups=5，表示有5组
warp="12"   # 
job_max=4


# 准备相图端点结构的OUTCAR
acnn_outcar2seed /path/to/OUTCAR /path/to/seed_name.res (ScZrB-ScB2-end-1.res)
```




<span style="font-size: 20px; color: lightblue;"> 3. 准备SEED文件

<span style="font-size: 20px; color: lightblue;"> 4. 提交任务

(1) 开始产生结构：RSS中创建文件夹（与RELAX/dyn_batch_relax中的src地址一致）
运行指令：
```shell        
airss.pl -build -max 5000 -seed ScZrB (即SEED_NAME)
```

(2) 提交结构预测到集群

