#!/bin/bash

# 获取当前目录路径
current_directory=$(pwd)



    # 查找匹配的 DISP 目录
    matching_disp_dirs=$(find "$current_directory" -type d -name "DIS.*")

    echo ${matching_disp_dirs}
    
    # 遍历匹配的 DISP 目录，并将文件拷贝到这些目录中
    for disp_dir in $matching_disp_dirs; do
        cd ${disp_dir}
        cp ../../INCAR ./
        cp ../../POTCAR ./
        cp ../../run.sh ./
        sbatch run.sh
        cd ../../
    done


