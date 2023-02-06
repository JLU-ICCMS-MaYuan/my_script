#!/bin/bash

fname=$1

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

IonSteps=`grep '!' $fname | wc -l`
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
Nat=`grep 'number of k points' -B 2 $fname |head -n 1|awk {'print($1)'}`
StruLine=`expr $Nat + 5 `