#!/bin/bash

# 检查是否提供了文件名
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# 读取文件名
file="$1"

# 计算BEGIN_CFG和END_CFG的数量
begin_count=$(grep -c "BEGIN_CFG" "$file")
end_count=$(grep -c "END_CFG" "$file")


if [ "$begin_count" -eq "$end_count" ]; then
    echo "The counts are equal: $begin_count"
else
    echo "The counts are not equal."
    echo "BEGIN_CFG count: $begin_count"
    echo "END_CFG count: $end_count"
fi
