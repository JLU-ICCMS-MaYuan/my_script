#!/bin/bash

outcar_file=$1
enthalpy=$(grep -a 'enthalpy' $outcar_file | tail -n 1 | awk '{print $ 5}')
# 读取文件并将第七行存储为变量
line=$(sed -n '7p' POSCAR)
# 将第七行按空格分割成数组
numbers=($line) # 使用一对圆括号()来指示这是一个数组。在圆括号内部，我们使用空格作为分隔符，将字符串分割成多个元素，然后将这些元素依次保存到数组中。# 例如，如果$line的值为1 2 3 4，那么执行numbers=($line)命令后，$numbers的值将变为一个包含四个元素的数组，即numbers[0]=1、numbers[1]=2、numbers[2]=3和numbers[3]=4。
# 初始化总数为0
total=0

# 遍历数组并将数字相加
for num in "${numbers[@]}"; do # ${numbers[@]}表示将数组numbers中的所有元素展开，并作为独立的参数传递给命令
# @符号表示展开所有元素，而$符号表示引用一个变量。 例如，如果numbers数组中有四个元素，即numbers[0]=1、numbers[1]=2、numbers[2]=3和numbers[3]=4，那么${numbers[@]}将展开为一个由空格分隔的字符串1 2 3 4。
  total=$(echo "$total + $num" | bc)
done

enthalpy_per_atom=$(echo "scale=8; ${enthalpy} / ${total}" | bc) #  设置scale参数，bc默认只保留整数部分，可能会导致结果不准确。

echo ${enthalpy_per_atom}  ${total}