#!/bin/bash

fname=$1
echo "Analyzed by QE post-processer code provided by Zicong"
echo "https://gitee.com/xczics/qe-post-processer"
echo "Checking output vc-relax output of $fname ..."
Error=`grep '%%%%%%%%%%%%%%' $fname`
if [ -n "$Error" ] # 如果Error对应的字符串变量长度不为0，则为真
then
        echo "Some Error happended in your task,please read output file yourself."
        exit 0
fi

echo "Checking if scf all be done..."
ScfStat=`grep -n 'convergence NOT achieved after' $fname` 
if [ -n "$ScfStat" ] # 如果ScfStat对应的字符串变量长度不为0，则为真
then
        echo "scf not converged"
        echo "$ScfStat"
        exit 0
fi

IonSteps=`grep '!' $fname | wc -l`  # 统计离子步数
IonStepConverged=`grep 'End final coordinates' $fname`
if [ -z "$IonStepConverged" ]
then
        echo "$IonSteps ion steps haved been calculated, BUT NOT CONVERGED"
        CONVERGED="False"
else
        echo "$IonSteps ion steps haved been calculated."
        echo "$IonStepConverged"
        CONVERGED="True"
fi

echo "Out put the last crystal structure..."
Nat=`grep 'number of k points' -B 2 $fname |head -n 1|awk {'print($1)'}`  # 获得胞内的原子数
#  grep -B 2 除了显示符合样式的那一行之外，并显示该行之前2行的内容
StruLine=`expr $Nat + 5 `
if [ -z "$2" ] # 如果string长度为0则为真
then
        grep 'CELL_' -A $StruLine $fname |tail -n `expr $StruLine + 1`
        # 除了显示符合范本样式的那一行之外，并显示该行之后的内容
else
        grep 'CELL_' -A $StruLine $fname |tail -n `expr $StruLine + 1` > $2 echo "The lastest structure in the file \"$2\""
        echo "INFO-FOR-Script: IonStep ( $IonSteps ), Nat ( $Nat ), IonConverged ( $CONVERGED ).awk should print \$4,\$8,\$12"
fi