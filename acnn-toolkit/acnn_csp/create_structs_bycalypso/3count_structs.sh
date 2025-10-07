#!/bin/bash
#
# 脚本功能: (最终版)
# 1. 遍历所有子目录。
# 2. 统计 'POSCAR_*' 文件数量。
# 3. 将未完成的目录路径写入 unfinished.dat。
# 4. (新) 将所有已完成的目录路径写入 finished.dat。
# 5. 在屏幕上输出详细的检查过程和最终统计。

# --- 初始化 ---
UNFINISHED_FILE="unfinished.dat"
FINISHED_FILE="finished.dat"
# 创建一个临时文件来存放已完成的目录列表
FINISHED_LIST_TMP="finished_list.tmp"

# 清空或创建输出和临时文件
> "$UNFINISHED_FILE"
> "$FINISHED_FILE"
> "$FINISHED_LIST_TMP"

echo "======== 开始扫描子目录结构完整性 (兼容模式) ========"
echo

# --- 主循环 ---
# 使用传统的 find | while 管道，这是最通用的写法
find . -mindepth 1 -maxdepth 1 -type d | sort -V | while IFS= read -r dir; do
    # 打印正在处理的目录名
    echo "--> 正在检查目录: $dir"

    # 统计该子目录中 `POSCAR_*` 文件的数量
    count=$(find "$dir" -maxdepth 1 -type f -name 'POSCAR_*' | wc -l)

    # 判断数量是否达到5000
    if [ "$count" -ge 5000 ]; then
        # 数量达标，将目录名写入临时文件
        echo "$dir" >> "$FINISHED_LIST_TMP"
        echo "    状态: 已完成 (找到 $count 个结构)"
    else
        # 数量不达标，将路径写入最终的输出文件
        echo "$dir" >> "$UNFINISHED_FILE"
        echo "    状态: 未完成 (找到 $count 个结构，已记录到 $UNFINISHED_FILE)"
    fi
    echo # 输出一个空行
done

# --- 循环结束后，进行最终统计和文件处理 ---

# 从文件中计算已完成和未完成的目录数量
finished_count=$(wc -l < "$FINISHED_LIST_TMP")
unfinished_count=$(wc -l < "$UNFINISHED_FILE")

# ★★★★★ 核心修改 ★★★★★
# 将记录了所有已完成目录的临时文件，直接重命名为最终的输出文件 finished.dat
mv "$FINISHED_LIST_TMP" "$FINISHED_FILE"


# --- 输出最终结果 ---

# 在屏幕上打印总结信息
echo "=========================================="
echo "检查完成！"
echo
echo "---------- 最终统计结果 ----------"
echo "已完成 (结构数 >= 5000) 的子目录数量: $finished_count"
echo "未完成 (结构数 < 5000) 的子目录数量: $unfinished_count"
echo "------------------------------------"
echo
echo "详情请查看以下文件:"
echo "  -> ${FINISHED_FILE} (列出了所有已完成目录的路径)"
echo "  -> ${UNFINISHED_FILE} (列出了所有未完成目录的路径)"
echo "=========================================="
