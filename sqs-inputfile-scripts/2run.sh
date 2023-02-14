#!/bin/bash

# 在运行完2run.sh命令后，会得到一个叫sqscell.out的文件, 这个文件里面存储的是包含352个原子的超胞的全部扩胞方式, 我们只需要取其中的一种方案。
# 清空sqscell.out里面的内容，然后重新定义结构形状让atat去定向寻找。
echo "" > sqscell.out

# 在sqscell.out里面指定某种扩胞方式，如果只指定一种扩胞方式，那么就在第一行写1. 然后在后面三行写好自定义的扩胞方式即可。

content="1
2 0 0
0 2 0
0 0 2"
echo "$content" > sqscell.out
