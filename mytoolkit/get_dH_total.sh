#!/bin/bash


enthalpy=$(grep enthalpy OUTCAR | tail -n 1 | awk '{print $ 5}')
# 读取文件并将第七行存储为变量

echo ${enthalpy} 