#!/bin/bash

# 清空日志文件
> fp_none.log
> fp_succeed.log
> fp_failed.log

# 遍历当前目录下的所有task.*目录
for dir in task.*; do
    if [ -d "$dir" ]; then
        # 检查是否有OUTCAR文件
        if [ ! -f "$dir/OUTCAR" ]; then
            echo "$dir" >> fp_none.log
            continue
        fi

        # 获取OUTCAR文件路径
        OUTCAR="$dir/OUTCAR"

        # 检查是否能grep出指定的内容
        grep -s -a "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" "$OUTCAR" >/dev/null
        free_energie=$?
        grep -s -a "General timing and accounting informations for this job" "$OUTCAR" >/dev/null
        timing_info=$?
        grep -s -a "in kB  \*" "$OUTCAR" >/dev/null
        in_kb=$?

        # 获取NELM步数和Iteration步数，确保是数字
        nelm_steps=$(grep -s -a "NELM   =" "$OUTCAR" | awk '{print $3}' | tr -d ';')
        iteration_steps=$(grep -s -a "Iteration" "$OUTCAR" | wc -l)

        # echo "$nelm_steps" "$iteration_steps"

        # 满足条件的认为是收敛
        if [ "$free_energie" -eq 0 ] && [ "$timing_info" -eq 0 ] && [ "$in_kb" -ne 0 ] && (( $(echo "$nelm_steps > $iteration_steps" | bc -l) )); then
            echo "$dir" >> fp_succeed.log
        else
            echo "$dir" >> fp_failed.log
        fi
    fi
done
