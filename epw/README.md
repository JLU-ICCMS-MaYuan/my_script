# EPW使用教程

## <span style="color:red">  EPW计算超导流程

### <span style="coler:lightgreen"> 0.提交任务的命令
```shell
# qe计算自洽
qe_main.py -i POSCAR -j bash -p 200 -w scf  scf -m mode=scf ecutwfc=80 ecutrho=960 kpoints_sparse='8 8 8'  degauss=0.02 execmd='mpirun -np 8' npool=4 queue=lhy nbnd=50

# qe计算能带
qe_main.py -i POSCAR -j bash -w band eletron -m mode=eleband  ecutwfc=80 ecutrho=960 execmd='mpirun -np 8' npool=4 queue=local kinserted=200 charge_density_dat='scf/tmp/H3S1.save/charge-density.dat' data_file_schema_xml='scf/tmp/H3S1.save/data-file-schema.xml' nbnd=50

# qe计算nscf
qe_main.py -i POSCAR -j bash -w epw  scf  -m mode=nscf ecutwfc=80 ecutrho=960 kpoints_dense='8 8 8'  degauss=0.02 execmd='mpirun -np 8' npool=4  charge_density_dat='scf/tmp/H3S1.save/charge-density.dat' data_file_schema_xml='scf/tmp/H3S1.save/data-file-schema.xml'  k_automatic=False wan=False  occupations=smearing  queue=lhy

# epw计算能带，其实是epw调用wannier计算能带, 要将计算的wannier的能带和qe计算的能带做对比
epw_main.py -i POSCAR -j bash -w epw epw_run -m mode=epw_eband execmd='mpirun -np 8' npool=8 dvscf_dir='./save' nbndsub=7 bands_skipped='1:12' dis_froz_min=-7.635175  dis_froz_max=23 dis_win_max=60  proj='H:s S:s S:p' queue=local nk='8 8 8' bands_skipped=1:1
# proj='Li:s;p Hf:s;p;d H:s'


# epw计算声子, 要将计算的声子freq.dat和qe计算的声子做对比
epw_main.py -i POSCAR -j bash -w epw epw_run -m mode=epw_phono  execmd='mpirun -np 8' npool=8 dvscf_dir='./save' nbndsub=7 bands_skipped='1:12' dis_froz_min=-7.635175  dis_froz_max=23 dis_win_max=60  proj='H:s S:s S:p'  nk='8 8 8' nq='4 4 4'  queue=lhy

# epw计算超导性质
epw_main.py -i POSCAR -j bash -w epw epw_run -m mode=epw_elph  execmd='mpirun -np 8' npool=8 dvscf_dir='./save' nbndsub=7 bands_skipped='1:12' dis_froz_min=-7.635175  dis_froz_max=23 dis_win_max=60  proj='H:s S:s S:p'  nk='8 8 8' nq='4 4 4' nkf='12 12 12' nqf='6 6 6' fsthick=0.4 degaussw=0.1 degaussq=0.5  queue=local
```

### <span style="color:lightgreen"> 1.自洽计算
```shell
mkdir 0.phonon
cd 0.phonon
mpirun -np N $QEBIN/pw.x -npool N < scf.in > scf.out
```

### <span style="color:lightgreen"> 2.声子计算
```shell
cd 0.phonon
mpirun -np N $QEBIN/ph.x  -npool N < ph.in > ph.out
```

### <span style="color:lightgreen"> 3.QE能带计算，其实也是一种非自洽计算，只不过k点是高对称路径的
```shell
mkdir -p tmp/Nb4H14.save
cp ../phonon/tmp/Nb4H14.save/charge-density.dat     tmp/Nb4H14.save/
cp ../phonon/tmp/Nb4H14.save/data-file-schema.xml   tmp/Nb4H14.save/
mpirun -np 48 /work/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <eleband.in> eleband.out
```

### <span style="color:lightgreen"> 4.准备QE动力学矩阵和受力
```shell
python pp.py
```

### <span style="color:lightgreen"> 5.非自洽计算, k点是均匀网格点
```shell
cp ../phonon/tmp/Nb4H14.save/charge-density.dat     tmp/Nb4H14.save/
cp ../phonon/tmp/Nb4H14.save/data-file-schema.xml   tmp/Nb4H14.save/
mpirun -np 48 /work/home/mayuan/software/qe-7.1/bin/pw.x  -npool 4 <nscf.in> nscf.out
```

### <span style="color:lightgreen"> 6. EPW计算声子
```shell
mpirun -np N $EPWBIN/epw.x -npool N < epw.in > epw.out
```

计算超导时相关的epw.in参数设置总览
```fortran
&inputepw
  ...
  ...
/
```

关于结构信息
```fortran
  prefix      = 'MgB2',
  amass(1)    = 24.305,
  amass(2)    = 10.811
```

关于文件读写
```fortran
  outdir      = 'tmp'   ! qe生成的tmp文件的路径
  dvscf_dir   = 'save'  ! pp.py保存的动力学矩阵和自洽势文件的路径
  
  wannierize  = .true.  ! 使用W90库调用计算万尼尔函数，并将旋转矩阵写入文件filukk文件。如果为False，旋转矩阵信息不自动产生，从filukk文件中读取。

  ep_coupling = .true.  ! 计算电声耦合
  elph        = .true.  ! 计算电声耦合系数
  ephwrite    = .true. ! 如果设置了Eliasberg = .true.，那么设置 ephwrite = .true.
                      !  ephwrite = .true. 会产生4个文件，这4个文件都需要在求解 Eliashberg equations 计算超导温度时被用到
                      ! ‘ephmatXX’  (XX: pool dependent files) files with e-ph matrix elements within the Fermi window (fsthick) on fine k and q meshes on the disk
                      ! ‘freq’ file contains the phonon frequencies
                      ! ‘egnv’ file contains the eigenvalues within the Fermi window
                      ! ‘ikmap’ file contains the index of the k-point on the irreducible grid within the Fermi window.

  epbwrite    = .true.  ! 粗糙布洛赫表示的电声耦合矩阵元和相关的数据（比如动力学矩阵）写入磁盘，存储为prefix.epb
  epbread     = .false. ! 粗糙布洛赫表示的电声耦合矩阵元和相关的数据从prefix.epb文件中读取
  
  epwwrite = .true.     ! 粗糙瓦尼尔表示的电声耦合矩阵元和相关的数据（比如动力学矩阵）写入磁盘，存储为epwdata.fmt和XX.epmatwpX
  epwread  = .false.    ! 粗糙瓦尼尔表示的电声耦合矩阵元和相关的数据从epwdata.fmt和XX.epmatwpX中读取。这个开关一般用于续算，
                        ! 要设置epwread = .true., 需要同时设置 kmaps = .true., 这是因为设置epwread = .true.时, kmaps相关信息不会重新计算产生，必须读取，所以必须设置kmaps = .true.。并且确保事先已经设置好 epwwrite = .true.
                        ! kmaps = .true. 用于从 prefix.kmap and prefix.kgmap 中读取k+q->k的散射情况。
  
  etf_mem     =  1      ! etf_mem = 1 所有致密的Bloch-space的电声耦合矩阵元都存在内存中, 这种方式更快, 此时IO更慢但是要求更少的内存
                        ! etf_mem = 2 在mode上对致密网格插补部分做了一个附加回路。通过设置“nmodes”可以进一步降低内存需求。
                        ! etf_mem = 3 多用于输运性质计算
                        ! etf_mem = 0 更少的io但是要求更多的内存
```

关于Wannier90
```fortran
  nbndsub     =  5,     ! 要使用的wannier函数的数量。
  num_iter    = 500     ! 为了最小化而传递给Wannier90的迭代次数
  dis_froz_max= 8.8     ! Wannier90的冻结状态
  proj(1)     = 'B:pz'               ! 在万尼90计算中使用的初始投影，比较简单的方法是：proj(1) = 'random'
  proj(2)     = 'f=0.5,1.0,0.5:s'
  proj(3)     = 'f=0.0,0.5,0.5:s'
  proj(4)     = 'f=0.5,0.5,0.5:s'
```

关于电声耦合计算
```fortran
  iverbosity  = 2         ! 2 = 罗嗦的输出针对超导部分
                          ! 3 = 罗嗦的输出针对电声耦合部分
  fsthick     = 0.4       ! 考虑自能δ函数的费米表面窗口宽度，单位是eV。缩小这个值可以减少自能量计算中包含的带数量。
  degaussw    = 0.10      ! Smearing in the energy-conserving delta functions in [eV]，默认0.025
  nsmear      = 1         ! 用于计算声子自能的不同smearings。
  delta_smear = 0.04      ! 声子自能计算时每一次额外的smearing的改变量，单位时eV

```

关于超导计算
```fortran
  ! eps_acustic = 2.0    ! 在进行电声耦合计算和a2f计算时的声子频率的下边界，单位时cm-1
  ! degaussq     = 0.5      ! 对q的所有电声耦合求和时的Smearing，起始值是0.5 ，单位是meV
  ! delta_qsmear = 0.05     ! 每次额外qsmearing的能量变化为0.05 meV 

  ! degaussq和delta_qsmear两个参数在 preifx.a2f 中最后几行的有详细写到，可以对应着检查。
  ! preifx.a2f中第2-11列是smearing=0.5, 0.55,...., 0.95的a2f， 
  ! preifx.a2f中第12-21列是smearing=0.5, 0.55,...., 0.95的2*a2f/w的积分

  nqstep       = 500      ! 用于计算a2f的步数

  eliashberg  = .true.    ! 启动超导计算的相关开关 
                          ! Note: To reuse .ephmat, .freq, .egnv, .ikmap files obtained in a previous run, one needs to set ep_coupling =.false., elph =.false., and ephwrite =.false. in the input file.
  laniso = .true.         ! 在虚轴上求解各向异性 Eliashberg equations，ephmat, .freq, .egnv, .ikmap 文件都必须提取准备好，它们由 ephwrite =.true.产生
  limag = .true.          ! 在虚轴上求解 Eliashberg equations
  lpade = .true.          ! 用 Padé 近似法将虚轴Eliashberg方程延续到实轴。这个开关要搭配 limag =.true.一起用

  conv_thr_iaxis = 1.0d-4 ! 虚轴Eliashberg方程迭代解的收敛阈值。

  wscut = 1.0      ! 单位eV，在求解Elisashberg 方程时的频率上限，必须搭配limag = .true.使用。如果设置了nswi为不为0的值，wscut被忽略。

  nstemp   = 1     ! 用于超导、传输、惯性等的温度点 数目。
                   ! 如果nstemp为空，或者等于temps(:)中的条目数，则使用temps(:)中提供的温度。
                   ! 如果nstemp>2并且temps(:)中只给出两个温度，则生成一个均匀间隔的温度网格，点之间的步长由(temps(2) - temps(1)) / (nstemp-1)给出。这个网格包含nstemp点。
                   ! nstemp不能大于50。
  temps    = 15.00 ! 用于超导、输运、惯性等的温度值，以开尔文为单位。
                   ! 如果没有提供temp，则temps=300和nstemp =1。
                   ! 如果提供两个温度，其中temps(1)<temps(2)和nstemp >2，则将temps转换为具有nstemp点的等间距网格。 其中temps(1)和temps(2)分别为最小值和最大值，在这种情况下，点的间隔根据(temps(2) - temps(1)) / (nstemp-1)设置。
                   !    例如 nstemp = 5 
                   !         temps = 300 500
                   ! 其它设置情况中，temps将被视为一个列表，直接使用给定的温度[temps = 17 20 30]。以这种方式提供的温度不能超过50个温度点。
  nsiter   = 500   ! 求解实轴或虚轴Eliashberg方程时自洽循环的迭代次数。

  muc     = 0.16   ! 有效库仑势在Eliashberg方程中的应用。
```

关于均匀网格选择
```fortran
  nk1         = 6 ! 粗电子网格的尺寸，对应于outdir中的nscf的计算和wfs。
  nk2         = 6
  nk3         = 6

  nq1         = 3 ! 粗声子网格的尺寸，对应于nqs列表。
  nq2         = 3
  nq3         = 3

  mp_mesh_k = .true. ! 如果.true.，设置在irr. wedge中设置精细的电子网，否则整个BZ中使用均匀网格

  nkf1 = 36 ! 如果filkf没有给出，则使用nkf指定的精细电子网格
  nkf2 = 36
  nkf3 = 36

  nqf1 = 18 ! 如果filqf没有给出，则使用nqf1指定的精细声子网格
  nqf2 = 18
  nqf3 = 18
 /
```

关于高对称路径的选择以及绘制EPW计算出的声子谱和能带, 其中path.dat可以通过get_hspp.py获得
```fortran
  band_plot = .true.
  filkf = 'path.dat'
  filqf = 'path.dat'
```
**特别注意:nkf1, nkf2, nkf3和filkf是互斥的一组参数, nqf1, nqf2, nqf3和filqf是互斥的一组参数**
## <span style="color:red"> 输出文件的说明

### <span style="color:lightgreen"> prefix.a2f
```shell
 w[meV] a2f and integrated 2*a2f/w for   10 smearing values
   # 频率      0.5展宽a2f       0.55展宽a2f                0.95展宽a2f      0.55展宽2*a2f/w积分值    0.90展宽2*a2f/w积分值  0.95展宽2*a2f/w积分值
   0.2220954   0.0000000        0.0000000    ..........    0.0000000        0.0000000     ........     0.0000000             0.0000000   
   0.4441908   0.0000000        0.0000000    ..........    0.0000000        0.0000000     ........     0.0000000             0.0000000  
   0.6662862   0.0000000        0.0000000    ..........    0.0000000        0.0000000     ........     0.0000000             0.0000000  
   ....
   ....
  110.8256080  0.0000000        0.0000000    ..........    0.0000000        0.0000000     .........    0.8100505             0.8100711
  111.0477034  0.0000000        0.0000000    ..........    0.0000000        0.0000000     .........    0.8100505             0.8100711
 Integrated el-ph coupling
  #            0.8099266   0.8099382   0.8099509   0.8099647   0.8099796   0.8099956   0.8100128   0.8100311   0.8100505   0.8100711
 Phonon smearing (meV)
  #            0.5000000   0.5500000   0.6000000   0.6500000   0.7000000   0.7500000   0.8000000   0.8500000   0.9000000   0.9500000
Electron smearing (eV)   0.1000000
Fermi window (eV)   0.4000000
Summed el-ph coupling    0.8098717

```

### <span style="color:lightgreen"> prefix.imag_aniso_XX
```shell
   #        w [eV]         Enk-Ef [eV]            znorm(w)       delta(w) [eV]           nznorm(w)
   #      频率              K-S本征态          准粒子重整化          超导能隙        正常态的准粒子重整化
    4.0608226305E-03    3.4713453311E-01    2.6773968092E+00    1.2423286333E-02    2.7824857836E+00
    4.0608226305E-03    3.4713466020E-01    2.6770752113E+00    1.2425568574E-02    2.7822774330E+00
    4.0608226305E-03   -1.8081323365E-01    1.5270752836E+00    2.9489459883E-03    1.5372999634E+00
    4.0608226305E-03    3.1137803949E-01    1.4372251813E+00    2.0242356203E-03    1.4455655124E+00
    4.0608226305E-03    3.5558684908E-01    2.6787011483E+00    1.2431313953E-02    2.7832553614E+00
    4.0608226305E-03    3.5558697618E-01    2.6782842340E+00    1.2433582157E-02    2.7828838557E+00
    4.0608226305E-03   -2.8508947876E-01    2.3777963007E+00    1.1155461384E-02    2.4547922278E+00
    4.0608226305E-03    3.6912900928E-02    2.4017656120E+00    1.1310260322E-02    2.4863893037E+00
    ..............
    ..............
```

### <span style="color:lightgreen"> prefix.pade_aniso_XX

该文件包含与MgB2中类似的信息。imag_aniso_XX，但超导间隙沿实轴，通过Pade近似得到。如果输入文件中的iverbosity=2，则写入该文件

```shell
    #        w [eV]         Enk-Ef [eV]        Re[znorm(w)]        Im[znorm(w)]   Re[delta(w)] [eV]   Im[delta(w)] [eV]
    2.2209540681E-04    3.4713453311E-01    2.6810423978E+00    1.6483256084E-05    1.2451399443E-02   -7.7264781883E-08
    4.4419081363E-04    3.4713453311E-01    2.6810720931E+00    3.2907646420E-05    1.2451667085E-02   -1.5418170907E-07
    6.6628622044E-04    3.4713453311E-01    2.6811215982E+00    4.9214380887E-05    1.2452113164E-02   -2.3040340743E-07
    8.8838162725E-04    3.4713453311E-01    2.6811909325E+00    6.5344821163E-05    1.2452737694E-02   -3.0558346696E-07
    1.1104770341E-03    3.4713453311E-01    2.6812801236E+00    8.1240557953E-05    1.2453540692E-02   -3.7937693122E-07

```

### <span style="color:lightgreen"> prefix.imag_aniso_gap0_XX
该文件包含各向异性超导隙的分布Δnk(ω=0) (eV)，由虚轴计算得到。


### <span style="color:lightgreen"> prefix.pade_aniso_gap0_xx
该文件包含与prefix.imag_aniso_gap0_XX相似的信息，但是这些超导能隙是通过Pade近似在实轴上计算得到的。

### <span style="color:lightgreen"> prefix.qdos_XX
这个文件包括了超导态准粒子的态密度，


## <span style="color:red"> 注意事项以及细节要求：

### 续算注意事项
```shell
Note: To reuse .ephmat, .freq, .egnv, .ikmap files obtained in a previous run, 
one needs to set ep_coupling =.false., elph =.false., and ephwrite =.false. in the input file.
```


### 并行方式
在使用EPW计算时，因为目前的版本还不支持平面波plancewaves G并行，所以必须设置并行核数等于并行k点数，即-np和-npool必须一致：
```shell
mpirun -np N $QEBIN/pw.x   -npool N < scf.in > scf.out
mpirun -np N $QEBIN/pw.x   -npool N < nscf.in > nscf.out
mpirun -np N $QEBIN/ph.x   -npool N < ph.in > ph.out
mpirun -np N $EPWBIN/epw.x -npool N < epw.in > epw.out
```

### Wannier90均匀撒点的方式

`/home/h240012/soft/qe-7.0/W90/utility/kmesh.pl`可以产生均匀的k点，使用方法：
```shell
perl kmesh.pl 4 4 4
```
上述脚本产生的k点是分数坐标的，需要再成倒格矢变为支教坐标

特别注意：六角格子的第一布里渊区是正六边形，但是kmesh.pl撒点时会在一个平行四边形里面撒点，这就会导致绘制电声耦合强度图不在第一布里渊区内。

### EPW也可以处理出各向同性的超导温度

```shell
fila2f = '' (Default)
Description: Input file with isotropic Eliashberg spectral function. 
The file contains the Eliashberg spectral function as a function of frequency in [meV]. 
This file can only be used to calculate the isotropic Eliashberg equations. 
In this case *.ephmat, *.freq, *.egnv, and *.ikmap files are not required.
''
```

## <span style="color:red"> 报错集锦
### <span style="color:lightgreen"> 1. Error in routine davcio (116): error writing file "../tmp/Nb4H14.epmatwe1"

原因是多个epw的计算任务公用同一个tmp文件中的波函数。
因为epw在计算的时候会在`outdir`指定的目录中创建一个目录叫`Nb4H14.ephmat`
这个目录中保存了`egnv`, `ephmat*`, `freq`, `ikmap`。
这几个文件都是都是`ephwrite=.true.`产生的。所以你如果公用的话，很有可能导致这个文件被不同的epw.in参数覆盖，导致文件内容混乱。

### <span style="color:lightgreen"> 2. Error in routine read_ephmat (1): nnk should be equal to nkfs

这个错误发生在`ephwrite=.false.`时，读取`$outdir\Nb4H14.ephmat\ephmat*`时发生的错误.

### <span style="color:lightgreen"> 3. coarse k-mesh needs to be strictly positive in 1st BZ
这是因为你的`prefix.save`中存储的波函数文件被覆盖了, 并不是epw.in要求的均匀网格nscf下的波函数文件. 

1. 仔细检查你提交任务的脚本是不是里面重复覆盖了波函数文件
2. 仔细检查你的nscf.in里面的k点的是不是用kmesh.pl生成的

### <span style="color:lightgreen"> 4. Error in routine elphon_shuffle_wrap (121): error opening crystal.fmt

这个错是卡在计算完所有的不可以q点之后，Bloch2wane之前。

我推测是因为epw.in没有完整跑完了， crystal.fmt应该是在epw.in没有完整跑完后生成的文件，    
所以如果用`epwwrite    = .false.,  epwread     = .true.`续算的话，需要读取这个文件又找不到只能报错。
所以说，不能续算。

### <span style="color:lightgreen"> 5. Error in routine epw_readin (1):  must use same w90 rotation matrix for entire run
这是因为你在一个完整的epw都没有算完的情况下就贸然开了`epwread = .true.`, 切记不要贸然打开。你都没有epwwrite，怎么epwread呢？

### <span style="color:lightgreen"> 6. Error in routine ephwann_shuffle (1): Error allocating epmatwe

<span style="color:lightyellow"> **这个错是卡在计算完所有的不可以q点之后，Bloch2wane之前。在修改参数开始跑之前，最好搞清楚怎么续算？
因为整个epw计算没有完成，停在了Bloch2wane之前，所以epbread和epwread不能同时打开, 如果倔强的你同时打开了，就会报错`Error in routine epw_readin (1): epbread cannot be used with epwread`。

只能按照如下方式续算**
```shell
elph        = .true.            no change
epbwrite    = .true.    --->    epbwrite    = .false.  
epbread     = .false.   --->    epbread     = .true.   
epwwrite    = .true.            no change
epwread     = .false.           no change
```

不清楚为什么报这个错，但是epw.out里面提到在输入文件里面加上`use_ws = .true.   ! Use the optimal Wigner-Seitz cell`可能会有帮助。
所以我加上了它，md，还是没跑过去。

下面是epw在`https://docs.epw-code.org/doc/GaN-II.html?highlight=use_ws`帖子里面提到了GaN-II计算epw时用到了use_ws这个参数。
```shell
Because our k and q coarse grid are so coarse, we need use_ws.
However it makes the calculation heavier (more vectors). 
If your fine grids are big enough (6x6x6 or more), you probably can use Gamma centered WS cell and use use_ws = .false.
```

通过检查EPW源代码发现了报错的位置：
```Fortran
!EPW/src/ephwann_shuffle.f90
...
...
    IF (etf_mem == 0) THEN
      ALLOCATE(epmatwe(nbndsub, nbndsub, nrr_k, nmodes, nqc), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating epmatwe', 1)
      ALLOCATE(epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_g), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating epmatwp', 1)
      epmatwe(:, :, :, :, :) = czero
      epmatwp(:, :, : ,: ,:) = czero
    ELSE
      ALLOCATE(epmatwe_mem(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating epmatwe_mem', 1)
      epmatwe_mem(:, :, :, :) = czero
    ENDIF
...
...
```

**那么这个时候问题有解了，只需要找到关于内存设置的地方就行。原来是因为我之前用的是etf_mem = 0，后来改成etf_mem = 1就行。**
```shell
# 下面是关于内存设置的一个重要开关：etf_mem

etf_mem = 0, then all the fine Bloch-space el-ph matrix elements are stored in memory (faster). 
etf_mem = 1, more IO (slower) but less memory is required. 
etf_mem = 2, an additional loop is done on mode for the fine grid interpolation part. This reduces the memory further by a factor “nmodes”. 
etf_mem = 3 is like etf_mem = 1 but the fine k-point grid is generated with points within the fsthick window using k-point symmetry (mp_mesh_k = .true. is needed) and the fine q-grid is generated on the fly. 
etf_mem = 3 is used for transport calculations with ultra dense fine momentum grids.
```

### <span style="color:lightgreen"> 7. qe报错：tetrahedra need automatic k-point grid

报这个错是因为在做qe非自洽计算时，你用了指定坐标的k点设置方式，但是这种设置方式不允许四面体方法计算非自洽。

## <span style="color:red"> 相关帖子
1. 关于电声耦合计算：https://xh125.github.io/2021/12/16/QE-epw/
2. 关于超导温度的计算：https://blog.csdn.net/lielie12138/article/details/127283037
3. EPW官方培训school的PDF文件：Tue.6.Lafuente.pdf(特别有帮助)，Wed.9.Mori(对于算超导特别有帮助) 

### <span style="color:lightgreen"> 8. Error in routine mem_size_eliashberg (1): Size of required memory exceeds max_memlt

通过检查`Size of allocated memory per pool: ~=    5.1757 Gb`可以知道每个pool分配到内存大小

然后设置`max_memlt = 8`可以增大每个pool分配的内存。

### <span style="color:lightgreen"> 9. Error: dis_spheres_first_wann is larger than num_bands-num_wann+1 Error: examine the output/error file for details

这是由于在计算wannier的时候设置的frozen窗口里面的能带数大于实即体系中的能带数，应该好好仔细检查你的nscf.in自洽是不是包含了足够多的nbnd。
比如：我的体系中CeSc2H24需要对41个wannier轨道进行投影，但是我在nscf.in中实际只用默认的nbnd计算了35个能带。需要重新计算nscf并设置nbnd=100

### <span style="color:lightgreen"> 10. Error in routine elphon_shuffle_wrap (1): Error allocating epmatq
这个是因为内存不够报错. 或者是设置nkf=1 1 1和nqf=1 1 1, 这导致k的电子态无法散射到k+q的电子态上插值出现问题。可以尝试增大fsthick

### <span style="color:lightgreen"> 11. Error in routine dynmat_asr (1):   wrong qpoint
这个是因为epw_phono.in里面设置的nq1, nq2, nq3的网格和qe设置的q网格不一致导致的。

### <span style="color:lightgreen"> 12. Error in routine loadqmesh_serial (1):  Cannot load fine q points
载入的nqf有问题

### <span style="color:lightgreen"> 13. coarse k-mesh needs to be strictly positive in 1st BZ
你可以仔细检查一下nscf.out中k点个数和epw_eband.out中的k点个数是不是一致，如果不一致，就修改nscf.in中的nscf为bands，重新非自洽。
这时候就可以得到一样数量的k点了。

### <span style="color:lightgreen"> 14. param_read_projection: malformed projection definition
这是是因为你用了,分割轨道信息，应该用分号
```shell
 proj(1)     = Li:s,p      --->   proj(1)     = Li:s;p  
 proj(2)     = Hf:s,p,d    --->   proj(2)     = Hf:s;p;d
 proj(3)     = H:s         --->   proj(3)     = H:s     
```

### <span style="color:lightgreen"> 15. Error in routine elphon_shuffle_wrap (1): Problem with modes file
