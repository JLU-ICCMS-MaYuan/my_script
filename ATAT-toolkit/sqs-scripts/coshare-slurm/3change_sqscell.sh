#!/bin/bash

# 清空sqscell.out里面的内容，然后重新定义结构形状让atat去定向寻找。
echo "" > sqscell.out


# 在sqscell.out里面指定你要定向寻找x种结构, 每种结构的超胞形式。

# Define the content
content="1
2 0 0
0 2 0
0 0 2"

# Write the content to the sqscell.out file
echo "$content" > sqscell.out
