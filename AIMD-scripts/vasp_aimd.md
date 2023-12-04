
# 1. 获得内能
```shell
grep "free  energy" OUTCAR|awk ' {print $5}' > energy.dat
```

# 2. 获得压强
```shell
grep "external pressure" OUTCAR | awk '{print $4}' > externalpressure.dat
```

在获得压强信息之后，我得解释一下vasp中OUTCAR的压强信息都代表什么含义？不然你个憨批啥也不知道！


## Pullay stress 是密度泛函理论（DFT）计算中的一种数值技术，用于修正晶格优化中的原子受力，数值上等于我们预设的压强（即参数PSTRESS就是Pullay stress的简称，巧妙吧，我也是刚刚发现P就是Pullay的意思），VASP默认为0（常压）。


```shell
# 这是分子动力学NVT中某一步的结果
FORCE on cell =-STRESS in cart. coord.  units (eV):
Direction    XX          YY          ZZ          XY          YZ          ZX     
--------------------------------------------------------------------------------------
Alpha Z  3740.33075  3740.33075  3740.33075
Ewald  -16397.91776-16401.71125-16355.31168   -63.42728   -46.44053   -55.84889
Hartree  2103.69118  2102.83931  2113.93792   -20.83328   -15.04682   -17.50973
E(xc)   -2383.87657 -2383.92845 -2383.68820    -0.12089    -0.15153    -0.13697
Local    3563.92027  3568.63607  3515.18307    77.72712    56.90523    67.42976
n-local  -565.15783  -564.85199  -566.92768    -1.27396    -0.84376    -0.66318
augment   652.63793   652.48765   653.02064     0.18755     0.06058     0.07825
Kinetic  9896.91377  9895.49527  9900.38867     2.80383     3.84077     2.39457
Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
-------------------------------------------------------------------------------------
Total     610.54174   609.29736   616.93349    -4.93693    -1.67605    -4.25620
in kB    1723.15215  1719.64010  1741.19180   -13.93365    -4.73037   -12.01240
external pressure =     1727.99 kB  Pullay stress =        0.00 kB                                                                                                                        

kinetic pressure (ideal gas correction) =     26.23 kB
total pressure  =   1754.22 kB
Total+kin.  1750.111    1745.427    1767.130     -14.300      -0.345     -13.106
```

在解释含义之前，你会发现一些关系：
external pressure = 1727.99 其实等于 1723.15215 + 1719.64010 + 1741.19180 = 1727.99
total pressure = external pressure + kinetic pressure =  1727.99 + 26.23 =  1754.22
Total in kB 在 XX, YY, ZZ, XY, YZ, ZX 6个方向上的分量就是一个晶胞感受到的应力，而这个应力的对角项 XX, YY, ZZ的和的平均值，就是我们施加在其上的外压力。


```shell
# 这是结构优化中某几步的结果
# 第一个离子步
  Total      42.96542    42.72859    43.06747     0.00000     0.00000     0.00000
  in kB     378.71198   376.62448   379.61151     0.00000     0.00000     0.00000
  external pressure =       78.32 kB  Pullay stress =      300.00 kB                                                            
  volume of cell :      181.77

......
# 第二个离子步
  Total      41.52877    41.33031    41.63164     0.00000     0.00000     0.00000
  in kB     363.70538   361.96725   364.60630     0.00000     0.00000     0.00000
  external pressure =       63.43 kB  Pullay stress =      300.00 kB                                                                                 
  volume of cell :      182.94

......
# 第三个离子步
  Total      37.37037    37.24460    37.42202     0.00000     0.00000     0.00000
  in kB     321.06699   319.98649   321.51073     0.00000     0.00000     0.00000
  external pressure =       20.85 kB  Pullay stress =      300.00 kB                                                                                 
  volume of cell :      186.48

```

在解释含义之前，你会发现一些关系：
Total in kB = 378.71198 + 376.62448 + 379.61151 = 378.31599
external pressure = Total in kB - Pullay stress = 78.32 

Total in kB 在 XX, YY, ZZ, XY, YZ, ZX 6个方向上的分量就是一个晶胞感受到的应力，而这个应力的对角相 XX, YY, ZZ的和就是我们施加在其上的外压力。

## external pressure在vasp做结构优化时的含义

external pressure如果大于0，则表示该结构在该压强下会膨胀，随着结构优化external pressure变小，晶格会随之变大，external pressure如果小于0，则表示该结构在该压强下会压缩，晶格会变小一些。那么这是为什么呢？？？

首先指明，external pressure其实是Total in kB 和 Pully Stress之间的差值，通过结构优化，external pressure最终会变为非常接近0的一个数，也就表明此时晶格感受到的外压和我们设置的外压基本一致。大白话理解是：通过结构优化将外压的效果转化为晶格参数的膨胀或者压缩。而衡量晶格此时所处的压强（即物理量Total in kB） 就是通过计算能量对晶格参数的二阶导获得的。

再来看上面的规律，
Total in kB - Pully Stress > 0, external pressure > 0, 也就意味着晶格感受到比你设置的压强更大的压强，所以需要释放原子之间的相互作用力，那么就需要膨胀晶格，增大原子间的距离，削弱原子间的受力。

Total in kB - Pully Stress < 0, external pressure < 0, 也就意味着晶格感受到比你设置的压强更小的压强，所以需要增加原子之间的相互作用力，那么就需要压缩晶格，增大原子间的距离，削弱原子间的受力。那么此时你一定会问一个问题，如果我把一个低压下的结构优化到高压下，我还能获得一个原子间受力为0的、处于势能面极小值点的结构吗？我想说，当然可以了！一个晶格受到的压强是通过计算能量关于晶格参数（其实就是一个量纲为m的物理量）的二阶导获得的，而原子间受力是通过计算能量关于原子坐标（也是一个量纲为m的物理量）的一阶导获得的，所以原子受力和晶格压强彼此没有任何关系，通过压缩晶格使之成为高压下的晶格并不影响晶格中的原子在结构优化后受力为零。

比如我下面放一个OUTCAR中grep出来的数据，你可以看到我在300kbar下做结构优化，一开始external pressure是一个比较大的数，但是当结构优化最终收敛时，它逐渐趋于零。

```shell
$ grep "external pressure" OUTCAR 
  external pressure =       78.32 kB  Pullay stress =      300.00 kB
  external pressure =       63.43 kB  Pullay stress =      300.00 kB
  external pressure =       20.85 kB  Pullay stress =      300.00 kB
  external pressure =        0.69 kB  Pullay stress =      300.00 kB
  external pressure =        0.13 kB  Pullay stress =      300.00 kB
  external pressure =        0.05 kB  Pullay stress =      300.00 kB
  external pressure =       -0.10 kB  Pullay stress =      300.00 kB
  external pressure =       -0.43 kB  Pullay stress =      300.00 kB
  external pressure =       -0.30 kB  Pullay stress =      300.00 kB
  external pressure =       -0.24 kB  Pullay stress =      300.00 kB
  external pressure =        0.01 kB  Pullay stress =      300.00 kB
  external pressure =        0.14 kB  Pullay stress =      300.00 kB
  external pressure =       -0.10 kB  Pullay stress =      300.00 kB
  external pressure =       -0.40 kB  Pullay stress =      300.00 kB
  external pressure =       -0.25 kB  Pullay stress =      300.00 kB
  external pressure =        0.19 kB  Pullay stress =      300.00 kB
  external pressure =        0.24 kB  Pullay stress =      300.00 kB
  external pressure =       -0.05 kB  Pullay stress =      300.00 kB
  external pressure =       -0.25 kB  Pullay stress =      300.00 kB
  external pressure =        0.00 kB  Pullay stress =      300.00 kB
  external pressure =        0.01 kB  Pullay stress =      300.00 kB
  external pressure =        0.09 kB  Pullay stress =      300.00 kB
  external pressure =       -0.02 kB  Pullay stress =      300.00 kB
  external pressure =       -0.04 kB  Pullay stress =      300.00 kB
  external pressure =       -0.02 kB  Pullay stress =      300.00 kB
  external pressure =       -0.01 kB  Pullay stress =      300.00 kB
```

## external pressure在做分子动力学时的含义

### external pressure在做NVT时的含义

首先明确在做NVT时，我们并不会设置PSTRESS，所以 Pullay stress = 0。那么在每次计算完一个离子步后输出的external pressure表示当前晶格感受到的外压是多少。按理说，如果你的结构是在100 GPa下优化好的，并且你这个结构在做NVT时保持稳定，那么计算出来的external pressure就会在100GPa附近上下波动。但是这里必须要指出：你并不能通过检查external pressure是否在某个压强点附近上下波动来判断你的结构是否稳定，在做NVT时，external pressure没有判断稳定性的价值。（这是我的好基友一个分子动力学大佬告诉我的，我也不知道为什么？）

# 3. 获得体积：
```shell
grep "volume of cell" OUTCAR | awk '{print $5}' > volume.dat
```

# 4. 获得温度：
```shell
grep "T= " OSZICAR | awk '{print $1 " " $3 " " }' > temperature_energy.dat
```


# 如何计算MSD?

## MSD的物理意义

单个粒子在单位时间间隔内的位移

## MSD的两种计算方法


## 什么是质心漂移？



# 如何计算RDF?

## MSD的物理意义

## RDF的计算方法