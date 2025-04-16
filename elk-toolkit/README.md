# ELK使用教程

## <span style="color:red"> 计算EOS并与其他第一性原理软件计算的EOS对比


## <span style="color:red"> 网络热帖子：
1. https://zhuanlan.zhihu.com/p/758415675
2. https://zhuanlan.zhihu.com/p/758415675
3. https://zhuanlan.zhihu.com/p/687280376#:~:text=ELK%20%E4%BB%A3%E7%A0%81%E4%B8%AD%E7%BC%96%E8%AF%91%E4%BA%86
4. https://sourceforge.net/p/elk/discussion/897820/thread/d5b59f97/#:~:text=I%20am%20a%20new%20user%20of%20Elk%20code.%20Now%20I
5. deepmd公司写的测试eos的：https://bohrium.dp.tech/notebooks/77351186918
6. ELK官网给出了一些非常重要的问题回答https://elk.sourceforge.io/faq.html。里面的第7个Basis set and convergence回答了很多重要的问题。
   1. 比如那些参数需要关注收敛性: rgkmax, gmaxvr, lmaxapw, lmaxvr, nempty, lradstp, swidth


## <span style="color:red"> 输入文件的参数解读

### 结构优化参数
```shell
latvopt = 0
# type of lattice vector optimisation to be performed during structural relaxation
# 在结构放松过程中晶格的形状和体积都将保持不变，系统只会优化原子的位置，而不会调整晶格参数。

latvopt = 1 
# the lattice vector optimisation will be constrained only by symmetry。
# 执行完全无约束的晶格向量优化。在这种模式下，晶格形状和体积, 原子位置都可以改变

latvopt = 2
# Optimisation over all symmetry-preserving strains except isotropic scaling is performed.
# It means an an iso-volumetric optimisation
# 进行等体积（iso-volumetric）优化. 这意味着在优化过程中虽然晶格形状, 原子位置都会改变，但总的体积不会改变。
```
## <span style="color:red"> 输出文件解读

### INFO.OUT中可以检查MPI运行是否正确

## <span style="color:red"> ELK自带小程序

### `spacegroup`可以产生晶体结构

### `eos`可以拟合equations of state to energy-volume data

### `findprimcell`

### `findsym`和`findsymcrys`

## <span style="color:red"> 报错集锦

### Warning(occupy): not enough empty states for k-point  

解决方法：increase nempty and also to perform a converge test of the total energy

### Warning(occupy): minimum eigenvalue less than minimum linearisation energy 
