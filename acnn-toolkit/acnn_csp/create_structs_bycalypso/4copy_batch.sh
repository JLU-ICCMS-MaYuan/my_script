#!/bin/bash

# --- 脚本设置 ---
# 定义要拷贝的源文件名
SOURCE_FILE="create_and_convert.sh"

# --- 1. 检查命令行参数 ---
# 检查用户是否提供了输入文件名作为参数
if [ -z "$1" ]; then
    echo "错误: 请提供一个包含目录列表的文件名作为参数。"
    echo "用法: $0 <文件名>"
    echo "例如: $0 finished.dat  或者  $0 unfinished.dat"
    exit 1
fi

# 将第一个命令行参数赋值给变量
INPUT_FILE="$1"

# --- 2. 检查所需文件是否存在 ---
# 检查源文件是否存在于当前目录
if [ ! -f "$SOURCE_FILE" ]; then
    echo "错误: 源文件 '$SOURCE_FILE' 在当前目录中不存在！"
    exit 1
fi

# 检查指定的输入文件是否存在且不为空
if [ ! -s "$INPUT_FILE" ]; then
    echo "错误: 指定的列表文件 '$INPUT_FILE' 不存在或为空！"
    exit 1
fi

# --- 3. 主循环 ---
echo "准备从列表文件 '$INPUT_FILE' 中读取目录并拷贝 '$SOURCE_FILE'..."
echo "--------------------------------------------------------"

# 使用 while read 循环，逐行读取指定的输入文件
# 这是处理文件内容最安全、最稳健的方法
while IFS= read -r dir_path || [ -n "$dir_path" ]; do
    # 检查从文件中读取的路径是否是一个真实存在的目录
    if [ -d "$dir_path" ]; then
        # 执行拷贝命令
        cp "$SOURCE_FILE" "$dir_path"
        # 打印一条信息，方便确认操作
        echo "已拷贝到 -> ${dir_path}"
    else
        # 如果不是目录，则打印警告信息并跳过
        echo "警告: '$dir_path' 不是一个有效目录，已跳过。"
    fi
done < "$INPUT_FILE"

echo "--------------------------------------------------------"
echo "拷贝完成！"
