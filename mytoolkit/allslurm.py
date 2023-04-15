#!/usr/bin/env python
import os

# 这是一个查看slurm任务运行情况的python脚本
#sq |grep R | awk 'NR!=1 {print $1, $5, $11}'
print("running")
tasks = os.popen('squeue -u liuhanyu -h --format "%A    %t    %R    %M    %Z     %s" | grep " R "').read().split('\n')
for t in tasks:
    if 'mayuan' in t:
        print(t)

print("waiting")
tasks = os.popen('squeue -u liuhanyu -h --format "%A    %t    %R    %M    %Z     %s" | grep " PD "').read().split('\n')
for t in tasks:
    if 'mayuan' in t:
        print(t)
