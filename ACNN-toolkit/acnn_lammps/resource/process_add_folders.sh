#!/bin/bash

# 查找所有以 _add 结尾的文件夹
for dir in *_add; do
  # 检查是否为目录
  if [ -d "$dir" ]; then
    echo "Entering directory $dir"
    cd "$dir" || exit

    # 执行指定的 Python 脚本
    python ~/playground/torchdemo/interface/calypso2/resource/convert/check_dt.py 100

    # 遍历当前目录中的所有文件
    for file in *; do
      # 检查是否是文件而不是目录
      if [ -f "$file" ]; then
        # 检查文件名是否包含 'xsf'
        if [[ "$file" != *xsf* ]]; then
          # 删除不包含 'xsf' 的文件
          echo "Deleting $file"
          rm "$file"
        fi
      fi
    done

    # 返回到上一级目录
    cd ..
  fi
done

