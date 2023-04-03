#!/bin/bash
read -p "指定开始的编号，例如输入1， 则从2开始编号" > i

for dir in *; do
    if [ -d "$dir" ]; then
        i=$((i+1))
        mv $dir ${i}
    fi
done

