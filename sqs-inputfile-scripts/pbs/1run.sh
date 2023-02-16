#!/bin/bash

# number of atoms = the number of atoms in primitive cell * dim = 44 * 8 = 352
 for id in {0..9}
 do
   mcsqs -n=352 -ip=$id &
   done
 wait


# 在运行完上面的命令后，会得到一个叫sqscell.out的文件。
# 当然也会产生bestcorr.out文件和bestsqs.out文件，但是我们主要目的还是获得sqscell.out文件。
# 不太清楚这个文件是干什么用的。但是一旦它出现以后，就将任务停掉。
