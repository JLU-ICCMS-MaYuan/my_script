# <div align="center"> <span style="color:red"> VASP 报错集合 </span> </div>

## VASP安装教程
https://zhuanlan.zhihu.com/p/547184583


## 自洽不收敛
### <span style="color:green"> 原子数特别多，元素种类特别多的自洽收敛  </span> </div>


达到步数上限时SCF仍没有收敛是实际研究中极为容易碰见的问题，VASP中经常会碰到电子不收敛的情况，也叫做SCF不收敛。经常会在电子结构比较复杂的体系中碰到。

1. 先检查几何结构是不是合理，非常离谱的初始结构会导致SCF难收敛。
2. 依次检查已经设置的参数，对照ppt的建议，是不是选择合理，对于+U和ICHARGE=11的任务，添加LMAXMIN=4 (对于d区体系)，添加LMAXMIN=6 (对于f区体系)。
3. 检查是不是ISTART=1读取了不合理的波函数，如果是，rm WAVECAR CHGCAR重新跑。
4. 如果ALGO=Fast或者VeryFast，换成Normal。
5. 尝试使用更大的SIGMA值，先粗略收敛，再读取CHGCAR和WAVECAR用小的SIGMA计算。
6. 对于非磁性体系（闭壳层ISPIN=1）添加：（注意AMIX和BMIX对收敛有很大影响，可以自己调试）
```shell
AMIX = 0.02
BMIX = 0.0001 # #almost zero, but 0 will crash some versions
```
7. 对于磁性体系（自旋极化，ISPIN=2）添加：
```shell
AMIX = 0.2
BMIX = 0.0001 #almost zero, but 0 will crash some versions
AMIX MAG = 0.8
BMIX MAG = 0.0001 #almost zero, but 0 will crash some versions
```
8. 尝试更换不同的ISMEAR
9. 检查体系是不是特殊的磁性排列，即MAGMOM设置是否合理
10. 提高积分精度，PREC=Accurate
11. 提高格点精度，ADDGRID = .TRUE. 
12. 先用1 1 1 K点计算收敛，再读取CHGCAR，用高K点计算
13. 尝试不同的ALGO 比如：ALGO=Conjugate
14. 如果在结构优化或者MD过程中，某一步突然不收敛，使用MAXMIX = 50
15. 尝试用更小的ENCUT或者更大的ENCUT的预收敛
16. 换更小的赝势或者更soft的赝势
17. 最后给出一个VASP官方教程里解决的不收敛的方法：
```shell
1. 用ALGO=N （是否收敛N to 2，Y to 6）
2. ICHARG=12 (no charge update，非自洽计算，N to5, Y to 3)
3. ICHARG=2 AMIX=0.1; BMIX=0.01(N to 4, Y to 6)
4. increase BMIXBMIX=3.0; AMIN=0.01(N to 5, Y to 6)
5. Bug report
```