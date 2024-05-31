#!/bin/bash
# 这是一个查看slurm任务运行情况的shell脚本
#sq |grep R | awk 'NR!=1 {print $1, $5, $11}'
echo "running"
squeue -u may -h --format '%A    %t    %R    %M    %Z     %s' | grep " R "
echo "waiting"
squeue -u may -h --format '%A    %t    %R    %M    %Z     %s' | grep " PD "



#a=`squeue | grep R | awk 'NR!=1 {print $1}'`
#for i in $a
#do echo $i
#scontrol show job $i | grep WorkDir 
#done

