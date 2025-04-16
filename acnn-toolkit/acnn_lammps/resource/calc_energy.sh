#!/bin/bash

# 检查是否传递了 "delet" 参数
delete_files=false
if [[ "$1" == "delet" ]]; then
    delete_files=true
fi

# 遍历当前路径下的所有 *_add 文件夹
for folder in *_add; do
    # 遍历文件夹中的所有 .xsf 文件
    for file in "$folder"/*.xsf; do
        if [[ -f "$file" ]]; then
            # 读取第一行中的 energy 值
            energy=$(grep "^# total energy" "$file" | awk '{print $5}')
            
            # 读取第11行中的原子数值 n
            n=$(sed -n '11p' "$file" | awk '{print $1}')
            
            # 计算 energy/n 的值
            result=$(echo "$energy / $n" | bc -l)
            
            # 判断结果是否大于 -6
            if (( $(echo "$result > -6.5" | bc -l) )); then
                echo "File: $file, Energy per atom: $result"
                
                # 如果传递了 "delet" 参数，删除文件
                if [[ "$delete_files" == true ]]; then
                    rm "$file"
                    echo "Deleted: $file"
                fi
            fi
        fi
    done
done

