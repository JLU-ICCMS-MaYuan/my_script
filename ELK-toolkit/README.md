# ELK使用教程

## <span style="color:red"> 计算EOS并与其他第一性原理软件计算的EOS对比


## <span style="color:red"> 网络热帖子：
https://zhuanlan.zhihu.com/p/758415675
https://zhuanlan.zhihu.com/p/687280376#:~:text=ELK%20%E4%BB%A3%E7%A0%81%E4%B8%AD%E7%BC%96%E8%AF%91%E4%BA%86
https://sourceforge.net/p/elk/discussion/897820/thread/d5b59f97/#:~:text=I%20am%20a%20new%20user%20of%20Elk%20code.%20Now%20I
deepmd公司写的测试eos的：https://bohrium.dp.tech/notebooks/77351186918


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

## `spacegroup`可以产生晶体结构

## `eos`可以拟合equations of state to energy-volume data

## `findprimcell`

## `findsym`和`findsymcrys`



