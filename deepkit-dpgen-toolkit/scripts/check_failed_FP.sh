#!/bin/bash

while read i; do
    element=$(sed -n '6p' $i/POSCAR)
    taskid=$(awk -F'/' '{print $NF}' <<< $i)  # 通过 <<< 传递变量给 awk
    echo $taskid $element
done < fp_failed.log
