# ACNN

##  <span style="font-size: 30px; color: lightgreen;"> 安装教程

可能需要的安装包准备： gcc(9以上版本)， ase，  pymatgen，  cmake (3以上版本)， qhull  openblas airss  xmgrace  libtorch torchdemo(ACNN) lammps(安装到torchdemo的interface里)

1. https://www.mtg.msm.cam.ac.uk/Codes/AIRSS;
2. https://plasma-gate.weizmann.ac.il/pub/grace/src/; 
3. https://www.openblas.net
 
###  <span style="font-size: 30px; color: red;">  编译流程：

####  <span style="font-size: 25px; color: blue;">  编译openblas
```shell
make -j8 USE_OPENMP=1
make install PREFIX=/...
```

如果在make的时候报错，比如缺少：

就写一个叫env_gcc.sh的脚本，内容包含：如何source你的gcc以及相应的lib库, 这里给出一个例子
```shell
source /public/env/gcc-9.2.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/gmp-6.1.2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mpfr-4.0.1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/mpc-1.1.0/lib
```

如果报错：`error while loading shared libraries: libmpfr.so.6: cannot open shared object file: No such file or directory`
要么你按照这个帖子自己安装一个gcc.9.2: https://www.zhihu.com/people/ma-yuan-94-83/posts. 要么你自己安装一个mpfr.

####  <span style="font-size: 25px; color: blue;">  安装airss

<span style="font-size: 20px; color: lightblue;"> 解压
```shell
$ tar -xzvf airss-v0.9.4.tgz

$ ls
bin  CONTRIBUTORS  examples  external  include  lib  LICENCE  makefile  README.md  src  test  VERSION.md
```
<span style="font-size: 20px; color: lightblue;"> 安装
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


####  <span style="font-size: 25px; color: blue;"> 安装grace

```shell
xmgrace: ./configure --prefix=... ; make -j8 ; make install
```

####  <span style="font-size: 25px; color: blue;"> 安装libtorch

PyTorch官网下载: https://pytorch.org/
    
![Download libtorch](picture/libtorch.jpg)
        
        解压后export即可
####  <span style="font-size: 25px; color: blue;"> 安装acnn

<span style="font-size: 20px; color: lightblue;"> 解压进入torchdemo

<span style="font-size: 20px; color: lightblue;"> 修改prefix.cmake

(除前五行需要修改，后面的全注释)，创建文件夹build，进入后cmake
```shell
cmake ..

make
```


####  <span style="font-size: 25px; color: blue;">  安装acnn和lmp_mpi的接口
<span style="font-size: 20px; color: lightblue;"> 进入torchdemo的interface/lammps

<span style="font-size: 20px; color: lightblue;"> 安装

自己下载个压缩包(lammps-2Aug2023.tar.gz)，更改build_lammps_interface.sh中的压缩文件名之后，激活intel进行编译（sh build_lammps_interface.sh build 核数），lmp_mpi在文件夹lammps-acnn/build里

#### 使用说明

1.  创建任务：deploy-x86-64 -p 30 -s ScZrB -b Sc-Sc=2.20,Sc-Zr=2.25,Sc-B=1.85,Zr-Zr=2.30,Zr-B=1.90,B-B=1.40 -n intel6240r_192

2.  修改参数：

    (1). AIRSS: dyn_gcs|dyn_grand修改原子体积VPA和FORMULA参数；

    (2). DFT: 准备POTCAR-元素 (用cat复制，不要用cp)，dyn_vasp_in修改INCAR参数,mkfdt修改res_file地址；

    (3). RELAX: dyn_batch_relax修改src地址和TASK参数；

    (4). SEED: 准备相图端点结构的OUTCAR，outcar2seed /path/to/OUTCAR /path/to/seed_name.res (ScZrB-ScB2-end-1.res)

3.  开始使用：

    (1) 训练集：AIRSS中创建文件夹（与RELAX/dyn_batch_relax中的src地址一致），airss输入文件如下：

![airss input](picture/airss_input.png)
        
        运行指令：airss.pl -build -max 5000 -seed ScZrB (即SEED_NAME)

#### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request


#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
