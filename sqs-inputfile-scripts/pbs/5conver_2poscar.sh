#!/bin/bash

# 需要自己安装一个格式转化器
mkdir -p stru/poscar
for i in {0..9}
do
    cd stru
    cp ../bestsqs$i.out  ./
    sqs2poscar bestsqs$i.out
    mv bestsqs$i.out-POSCAR poscar/POSCAR-$i
    sed -i 's/xxx/1/g' poscar/POSCAR-$i
    cd ..
done

