#!/bin/bash

# 初始化计数变量，从 6 开始
x=0

# 目标文件夹
target_folder="40_add"

# 遍历 63_add 到 68_add 文件夹中的 *.xsf 文件
for folder in {41..62}_add; do
    for file in "$folder"/*.xsf; do
        if [[ -f "$file" ]]; then
            # 生成新文件名并拷贝到目标文件夹
            new_file="$target_folder/40_labeling_add_fp${x}_OUTCAR_struct000000.xsf"
            cp "$file" "$new_file"
            echo "Copied $file to $new_file"
            
            # 增加计数
            ((x++))
        fi
    done
done

