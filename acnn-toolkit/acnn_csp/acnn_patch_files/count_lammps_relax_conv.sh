#!/bin/bash

# 检查是否传入了一个参数
if [ $# -ne 1 ]; then
    echo "Usage: $0 <integer>"
    exit 1
fi

# 确保参数是整数
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
    echo "Error: The argument must be an integer."
    exit 1
fi

count=0
IT="$1"

# 使用 find 命令来计数
count=$(find . -maxdepth 3 -type f -wholename "./IT${IT}/AlBe*/*-out.res" | wc -l)

echo "Found $count directories starting with 'IT$IT'."

