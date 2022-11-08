#!/bin/bash

# number of atoms = the number of atoms in primitive cell * dim = 44 * 8 = 352
 mcsqs -n 352

# 在运行完上面的命令后，会得到一个叫sqscell.out的文件。
# 不太清楚这个文件是干什么用的。但是一旦它出现以后，就将任务停掉。
