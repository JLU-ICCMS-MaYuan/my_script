#!/bin/bash

# number of atoms = the number of atoms in conventional cell * dim = 44 * 8 = 352
# 其实这里的352限制的是扩胞后的胞内原子个数。
# 这里一定是晶胞的原子数*扩胞方式。
 mcsqs -n 352

# 在运行完上面的命令后，会得到一个叫sqscell.out的文件, 这个文件里面存储的是包含352个原子的超胞的全部扩胞方式, 我们只需要取其中的一种方案。

# 不太清楚这个文件是干什么用的。但是一旦它出现以后，就将任务停掉。
