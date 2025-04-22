# Wannier90使用教程

wannier90按照Marzari和Vanderbilt (MV)的方法计算最大局部万尼尔函数maximally-localised Wannier functions (MLWF)

##  <span style="color:red">  Wannier90拟合能带的基本流程

##  <span style="color:red">  用pw.x运行计算自洽（scf）和非自洽（nscf）

自洽网格与非自洽网格可以不一致。但是非自洽网格与`Ce1Sc2H24.win`中网格必须一致。

非自洽计算时必须手动设置好以下参数：
```Fortran
&control
 verbosity = 'high' #输出本征值     
/
&system
 nbnd = xxx !设定能带数, 千万要注意，能带数一定要多一点. 比如CeSc2H24, 用41个num_wannier, 那么能带数一定要超过41
/
&electrons
 diago_full_acc = .true. #对角化
/
```
非自洽生成k点的方法：
```shell
kmesh.pl 8 8 8  > kpoint 
```

##  <span style="color:red">  准备Ce1Sc2H24.win输入文件

```shell
begin projections
 Ce:f
 Sc:d
 H:s
end projections
exclude_bands = 1:29
                                                                                
begin unit_cell_cart
Ang
    5.2292802927    0.0000000000    0.0000000000
   -2.6146401463    4.5286895771    0.0000000000
    0.0000000000    0.0000000000    3.6396561976                                
end unit_cell_cart

begin atoms_frac
Ce  0.000000000000  0.000000000000  0.000000000000
Sc  0.666666666570  0.333333333141  0.500000000288
Sc  0.333333333156  0.666666666513  0.500000000288
H   0.237917437763  0.762082561912  0.000000000000
H   0.524165124248  0.762082561912  0.000000000000
H   0.762082562043  0.524165124285  0.000000000000
H   0.237917437684  0.475834875369  0.000000000000
H   0.762082562063  0.237917437742  0.000000000000
H   0.475834875479  0.237917437742  0.000000000000
H   0.000000000000  0.617417740465  0.217157425327
H   0.617417740594  0.617417740465  0.217157425327
H   0.000000000000  0.382582259189  0.217157425327
H   0.382582259233  0.382582259189  0.217157425327
H   0.382582259277  0.000000000000  0.217157425327
H   0.617417740723  0.000000000000  0.217157425327
H   0.773679796493  0.773679796343  0.500000000288
H   0.226320203334  0.226320203310  0.500000000288
H   0.773679796642  0.000000000000  0.500000000288
H   0.226320203358  0.000000000000  0.500000000288
H   0.000000000000  0.226320203310  0.500000000288
H   0.000000000000  0.773679796343  0.500000000288
H   0.382582259233  0.382582259189  0.782842575249
H   0.000000000000  0.382582259189  0.782842575249
H   0.000000000000  0.617417740465  0.782842575249
H   0.382582259277  0.000000000000  0.782842575249
H   0.617417740594  0.617417740465  0.782842575249
H   0.617417740723  0.000000000000  0.782842575249
end atoms_frac
                    
num_wann =  41
iprint =   2
dis_froz_min =     8.000000000000
dis_froz_max =    32.620000000000
num_iter =     500
write_bvec = .true.

bands_plot = .true.
begin kpoint_path
G  +0.0000 +0.0000 +0.0000 M  +0.5000 +0.0000 +0.0000
M  +0.5000 +0.0000 +0.0000 K  +0.3333 +0.3333 +0.0000
K  +0.3333 +0.3333 +0.0000 G  +0.0000 +0.0000 +0.0000
G  +0.0000 +0.0000 +0.0000 A  +0.0000 +0.0000 +0.5000
A  +0.0000 +0.0000 +0.5000 L  +0.5000 +0.0000 +0.5000
L  +0.5000 +0.0000 +0.5000 H  +0.3333 +0.3333 +0.5000
H  +0.3333 +0.3333 +0.5000 A  +0.0000 +0.0000 +0.5000
end kpoint_path
bands_plot_format = gnuplot
num_bands = 50
dis_num_iter = 1000
write_hr = .true.

mp_grid = 6 6 6
begin kpoints
   0.00000000  0.00000000  0.00000000 
   0.50000000  0.00000000  0.00000000
   ...
   ...
End Kpoints

```


相关参数设置的经验贴：
1. wannier90拟合能带能量窗口参数说明：https://blog.csdn.net/bubu789/article/details/119220576
2. wannier90计算流程说明：https://zhuanlan.zhihu.com/p/381615718
3. 分享一个确定Wannier90能量窗口的脚本: https://blog.sciencenet.cn/blog-2909108-1263724.html
4. wannier90拟合能带能量窗口调节: https://zhuanlan.zhihu.com/p/541333688
5. wannier90参数说明：https://yxli8023.github.io/2021/08/03/Wannier90-Study.html
6. wannier90计算流程说明：https://yxli8023.github.io/2021/04/22/Wannier90-Band.html
7. https://github.com/Xin-Jiaqi/heterojunction/tree/main


**最重要的就是搞清楚如下参数的设置：**

#### 晶格，单位默认是angstrom
```shell
begin unit_cell_cart
Ang
     4.1137419799060595    0.0000000000000000    0.0000000000000000
     0.0000000000000000    4.1137419799060595    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.1137419799060595
end unit_cell_cart
```

#### 坐标，单位默认是angstrom
```shell
begin atoms_cart
Ang
Nb       0.7500000000000000  0.7500000000000000  0.7500000000000000
Nb       0.2500000000000000  0.2500000000000000  0.7500000000000000
Nb       0.2500000000000000  0.7500000000000000  0.2500000000000000
Nb       0.7500000000000000  0.2500000000000000  0.2500000000000000
H        0.6549958318805054  0.6549958318805054  0.3450041681194946
H        0.6549958318805054  0.3450041681194946  0.6549958318805054
H        0.3450041681194946  0.6549958318805054  0.6549958318805054
H        0.8450041681194946  0.8450041681194946  0.1549958318805054
H        0.1549958318805054  0.1549958318805054  0.1549958318805054
H        0.8450041681194946  0.1549958318805054  0.8450041681194946
H        0.1549958318805054  0.8450041681194946  0.8450041681194946
H        0.3450041681194946  0.3450041681194946  0.3450041681194946
H        0.5000000000000000  0.0000000000000000  0.0000000000000000
H        0.0000000000000000  0.5000000000000000  0.0000000000000000
H        0.0000000000000000  0.0000000000000000  0.5000000000000000
H        0.5000000000000000  0.0000000000000000  0.5000000000000000
H        0.0000000000000000  0.5000000000000000  0.5000000000000000
H        0.5000000000000000  0.5000000000000000  0.0000000000000000
end atoms_cart
```

#### KPOINTS 这是设置的网格必须与非自洽`nscf.in`中的网格一致
```
! KPOINTS
mp_grid : 8 8 8
```

#### 将每一个K点的坐标都写下来

可以利用qe安装目录下的`wannier90-3.1.0/utility/kmesh.pl`实现
```shell
# wan可以阻止输出权重
./kmesh.pl 8 8 8 wan
# 没有wan可以输出权重，用于写nscf.in中的ATOMIC_POSITIONS (crystal)
./kmesh.pl 8 8 8 
```

```shell
# 下面是`Ce1Sc2H24.win`中的kpints例子
begin kpoints
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.12500000
  0.00000000  0.00000000  0.25000000
  ....
  ....
  0.87500000  0.87500000  0.62500000
  0.87500000  0.87500000  0.75000000
  0.87500000  0.87500000  0.87500000
end kpoints
```

#### 解纠缠时设计的两个窗口`outer_window`和`frozen_window`

`outer_window`由`dis_win_min`, `dis_win_max`控制，在这两个值范围限制的能量范围内，至少包含num_wann条能带，最好费米能级在这个能量范围内。

实际上一般将outer_window设置为用exclude_bands排除无关能带后，剩下的能带能量范围，两端再留些余量。

`frozen_window`由`dis_froz_min`, `dis_froz_max`控制，在这两个值范围限制的能量范围内，设定为num_wann对应能带的能量范围，两端留些余量。

关于什么情况下需要解纠缠？简单来说：只有在拟合完整的半导体或绝缘体价带时，才不需要解纠缠。对于半导体的导带或金属能带，很难找到一组和其它能带完全分离的能带。这种情况下就需要解纠缠。具体看帖子：https://blog.csdn.net/bubu789/article/details/119220576

那么具体到底怎么设置`dis_win_min`, `dis_win_max`和`dis_froz_min`, `dis_froz_max`这四个值呢？

你可以同时打开`nscf.out`和`Ce1Sc2H24.win`。然后看着nscf.out中任意一个k点的能级，取大致估算你如何设置这四个窗口值。
```shell
          k = 0.8750 0.8750 0.5000 (  5670 PWs)   bands (ev):                  |
   -31.8058 -31.8046 -31.5983 -31.5972  -9.0731  -8.9430  -8.7745  -8.6940     |num_wann  = 34  
    -8.6485  -8.5790  -8.0634  -8.0351  -7.7716  -7.7565  -7.7516  -7.7414     |num_bands = 40
    11.0499  11.1362  13.0502  13.1029  14.4198  14.4694  14.6088  15.3367     |exclude_bands = 1-16
    15.7313  16.0706  17.5225  18.1113  19.4177  19.8881  22.1927  22.5352     |
    23.0776  23.1916  25.5394  25.5405  25.5776  26.3891  26.9640  27.7693     |
    28.9978  29.2331  29.3141  29.7651  30.8552  31.0280  31.4839  31.9324     |dis_win_min = 15 # eV
    32.0483  32.6754  33.5966  33.6699  36.4114  37.3127  37.5861  37.6481     |dis_win_max = 41 # eV
    39.1632  40.1968  40.3752  40.5535                                         |dis_froz_min = 15 # eV
                                                                               |dis_froz_max = 40 # eV
# 我们已经通过scffit.out知道Ce1Sc2H24的费米能级在23.4468eV.所有必须确保费米能级在这其中。
# 所以以22eV为中心，上下各取num_wann/2数量的能带即可，不用严格1/2, 稍微有点误差也行。
# 甚至你可以用60-34=26, 直接放弃最低能级的16条带。这很粗暴，一般情况奏效。
# 一定要确保num_bands + exclude_bands = qe中的nbnd。
# 检查Ce1Sc2H24.amn，确保其中的能带数等于Ce1Sc2H24.win中的num_bands
```

## <span style="color:red">  执行`wannier90.x -pp Ce1Sc2H24`获得Ce1Sc2H24.nnkp

执行这句话，wannier90.x会自动读取`Ce1Sc2H24.win`里面关于`num_wann`和`begin projections...end projections`的设置。

**所以如果你想修改`num_wann`和`begin projections...end projections`重新计算wannier90的能带，必须从这一步开始执行。**

`-pp`后面跟的是你自定义的prefix, 叫什么都可以，自己记着点就行。注意这个破程序`wannier90.x`不能并行。
```shell
wannier90.x -pp Ce1Sc2H24
```

### <span style="color:yellow">  执行` ~/software/qe-7.1/bin/pw2wannier90.x -pd .true. < Ce1Sc2H24.pw2wan.in > Ce1Sc2H24.pw2wan.out`获得`Ce1Sc2H24.mmn`, `Ce1Sc2H24.eig`, `Ce1Sc2H24.amn`.

`pw2win.inp`的文件模板
```shell
&inputpp
  outdir     =  './' # 如果设置'./', 程序就会在当前目录下寻找Ce1Sc2H24.save这个目录下的波函数文件。如果设置'../', 程序就会在上级目录下寻找Ce1Sc2H24.save这个目录下的波函数文件。
  prefix     =  'Ce1Sc2H24'       
  seedname   =  'Ce1Sc2H24'
  write_amn  =  .true.
  write_mmn  =  .true.
/
```



1. `Ce1Sc2H24.mmn` 重叠矩阵
2. `Ce1Sc2H24.amn` Bloch states 到一个局域轨道的投影
3.  `Ce1Sc2H24.eig` 每一个k点的bloch本征态

注意这个破程序`pw2wannier90.x`竟然可以并行。但是在wannier.90的手册中提到`Note that, unless you specify wf_collect=.true. in your pw.x input file, you must run pw2wannier90 with the same number of processors as pw.x`
```shell
pw2wannier90.x < pw2win.inp > pw2win.out
```

### <span style="color:yellow">  执行`wannier90.x Ce1Sc2H24.win > wannier90.log 2>&1`开始拟合wannier90能带

`wannier90.x 可以并行但是要在编译的时候搞好了`
```shell
# postw90.x and wannier90.x can be run in parallel to speed up the calculations, using the MPI libraries.

# To enable the parallel version to be built, you must specify some flags in the make.inc file of wannier90 and postw90; for further information, please refer to the README.install file in the top directory of the wannier90 distribution.

wannier90.x Ce1Sc2H24.win > wannier90.log 2>&1
```

##  <span style="color:red">  Wannier90计算报错集锦
### 1. 执行`wannier90.x -pp Ce1Sc2H24`报错`param_get_projection: Problem reading m state into string Error: examine the output/error file for details`

可能发生报错的位置：
1. `Ce1Sc2H24.win`里面`num_wann = 34`与`begin projections ... end projections`矛盾。比如下面的设置是合理的：
    ```shell
    num_wann  = 34  # 设置需要投影的Wannier轨道，4个Nb，每个Nb有5个d轨道，总共20个d轨道，14个H有14个s轨道，总共有34个轨道。所以需要设置num_wann轨道为34

    begin projections # 这是第1种写法
    Nb : dz2; dx2-y2; dxy; dyz; dxz 
    H  : s
    end projections

    begin projections # 这是第2种写法
    Nb : d 
    H  : s
    end projections
    ```
2. 千万注意`dz2; dx2-y2; dxy; dyz; dxz`的`dxz`后面不要加`;`。加了一定报错。

### 2.执行`wannier90.x Ce1Sc2H24.win > wannier90.log 2>&1`报错` dis_windows: More states in the frozen window than target WFs`

1. 报这个是因为你设置的`frozen_window`内的能带数超出了`num_wann`个数了。
2. 如果你是按照我前面讲的设置`frozen_window`的方法设置的`dis_froz_min`, `dis_froz_max`。那么其实你不需要大动干戈，只需要一个eV一个eV的减小`dis_froz_max`或者增加`dis_froz_min`即可缩小冻结在窗口中的轨道数。

### 3.执行`wannier90.x Ce1Sc2H24.win > wannier90.log 2>&1`报错`dis_windows: Energy window contains fewer states than number of target WFs`

1. 报这个是因为你设置的`dis_windows`内包含的态小于`num_wann`个数, 可以适当扩展`dis_win_min`和`dis_win_max`确定的能量范围。

### 4. 执行`wannier90.x Ce1Sc2H24.win > wannier90.log 2>&1`报错`too many projections to be used without selecting a subset`

### 5. 执行`wannier90.x -pp Ce1Sc2H24`  报错 `Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid Error: examine the output/error file for details`
这明显就是没有设置mp_grid

### 6. 执行`wannier90.x < Ce1Sc2H24.win > Ce1Sc2H24.wout` 报错 `param_read: mismatch in Ce1Sc2H24.eig Error: examine the output/error file for details`  
这明显就是