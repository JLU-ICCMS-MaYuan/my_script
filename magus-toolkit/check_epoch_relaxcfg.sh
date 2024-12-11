#!/bin/bash

# 输出文件
output_file="epoch_info.txt"
# 清空输出文件
> "$output_file"

# 写入表头
echo -e "Epoch Number    To_Relax_CFG_Count    Relaxed_CFG_Count" > "$output_file"

# 遍历所有以 epoch 开头的目录
for dir in epoch*/; do
  # 提取目录名中的epoch编号
  epoch_number=$(basename "$dir" | sed 's/epoch//')

  # 进入该目录
  cd "$dir" || continue

  # 检查文件是否存在并统计 END_CFG 出现次数，否则赋值为 Null
  if [ -f "to_relax.cfg" ]; then
    to_relax_count=$(grep -c 'END_CFG' to_relax.cfg)
  else
    to_relax_count="Null"
  fi

  if [ -f "relaxed.cfg" ]; then
    relaxed_count=$(grep -c 'END_CFG' relaxed.cfg)
  else
    relaxed_count="Null"
  fi

  # 将结果写入输出文件，使用大空格分隔列
  echo -e "$epoch_number              $to_relax_count                   $relaxed_count" >> "../$output_file"

  # 返回上一级目录
  cd ..
done

echo "结果已写入 $output_file"

