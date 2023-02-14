#!/bin/bash

# 需要自己安装一个格式转化器
# 安装好以后运行
# 这里稍微解释一下bestsqs.out文件
# 第一行到第三行是体系的晶格参数。第四行到第六行是体系的扩胞方案
sqs2poscar bestsqs.out 

# 得到一个名为bestsqs.out-POSCAR的文件
