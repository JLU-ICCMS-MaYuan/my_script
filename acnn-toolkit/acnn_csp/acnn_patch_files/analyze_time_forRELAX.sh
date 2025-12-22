#!/bin/bash

# 检查命令行参数是否传入
if [ $# -ne 1 ]; then
    echo "用法：$0 <包含Wall time的文件路径>"
    echo "示例：$0 IT31/relax_group_aa"
    exit 1
fi

# 从命令行参数获取文件路径
TARGET_FILE="$1"

# 检查文件是否存在
if [ ! -f "$TARGET_FILE" ]; then
    echo "错误：文件 '$TARGET_FILE' 不存在！"
    exit 1
fi

# 提取Wall time行（核心命令，适配参数传入的文件）
WALL_TIME_LINES=$(grep "Wall time" "$TARGET_FILE")
echo "$WALL_TIME_LINES"

# 检查是否提取到有效数据
if [ -z "$WALL_TIME_LINES" ]; then
    echo "错误：文件 '$TARGET_FILE' 中未找到包含'Wall time'的行！"
    exit 1
fi

# 初始化总秒数（浮点数）
TOTAL_SECONDS=0

# 遍历每行提取秒数并累加（用bc处理浮点数运算）
echo "$WALL_TIME_LINES" | while read -r _ _ sec _; do
    TOTAL_SECONDS=$(echo "$TOTAL_SECONDS + $sec" | bc -l)
done

# 转换为小时（1小时=3600秒）
TOTAL_HOURS=$(echo "$TOTAL_SECONDS / 3600" | bc -l)

# 格式化输出（保留3位小数，提升可读性）
FORMATTED_SECONDS=$(printf "%.3f" "$TOTAL_SECONDS")
FORMATTED_HOURS=$(printf "%.3f" "$TOTAL_HOURS")
LINE_COUNT=$(echo "$WALL_TIME_LINES" | wc -l)

# 输出最终结果
echo -e "\n===== Wall Time 统计结果 ======"
echo "目标文件：$TARGET_FILE"
echo "有效数据行数：$LINE_COUNT 行"
echo "总耗时（秒）：$FORMATTED_SECONDS s"
echo "总耗时（小时）：$FORMATTED_HOURS h"
echo "===============================\n"
