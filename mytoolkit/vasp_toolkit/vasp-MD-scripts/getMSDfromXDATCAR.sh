#!/bin/bash

# 输出提示信息
echo "NOTE --------------------"
echo "    Use vaspkit to get the MSD for each element provided"

# 检查是否有传入参数
if [ "$#" -eq 0 ]; then
  echo "No elements provided. Exiting..."
  exit 1
fi

# 遍历每一个传入的元素名称
for element in "$@"; do
  echo "Processing MSD calculation for element: $element"
  
  # 构建输入内容
  input="722\n 1\n $element\n 2000\n 1\n"
  
  # 打印输入内容并传递给 vaspkit
  echo -e "$input" | tee /dev/tty | vaspkit > /dev/null
  
  # 检查是否生成了 MSD.dat 文件，并重命名
  if [ -f "MSD.dat" ]; then
    mv "MSD.dat" "${element}_MSD.dat"
    echo "Renamed MSD.dat to ${element}_MSD.dat"
  else
    echo "MSD.dat not found for $element. Skipping rename."
  fi
  
  # 提示每个元素的处理完成
  echo "MSD calculation for $element completed"
done

echo "All MSD calculations finished."