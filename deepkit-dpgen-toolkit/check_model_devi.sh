#!/bin/bash  
 
beginid=$1  
endid=$2  

echo "You can use it by: ./check_model_devi.sh 000 820"
# 使用循环检查每个目录  
for i in $(seq -w "$beginid" "$endid"); do  
    directory="task.$i.000"  
    file="$directory/model_devi.out"  
 
    if [ -f "$file" ]; then  # 检查文件是否存在  
        count=$(grep -c 'step' "$file")  # 使用 -c 选项直接计算匹配行数  
        if [ "$count" -ne 1 ]; then  
            echo "$directory/model_devi.out contains $count 'step' lines."  
        fi  
    else  
        echo "File $file does not exist."  
    fi  
done