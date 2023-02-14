#!/bin/bash

mcsqs -rc

# 上面的命令运行后，着重看bestcorr.out bestsqs.out 两个文件。

#  bestsqs.out 是输出的结构文件
#  bestcorr.out 是判断结构好坏的判断文件。最后一列0越多，结构越好。objective_function的负的绝对值越大结构越好。
#  mcsqs.log 文件中的objective_function越小说明产生的结构越好

