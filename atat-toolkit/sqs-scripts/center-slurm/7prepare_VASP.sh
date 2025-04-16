#!/bin/sh 


# 获取当前工作目录
current_dir=$(pwd)

# 定义变量 potcar_lib
potcar_lib="potcar_lib"

# 构建新的路径
pot_PATH="${current_dir}/${potcar_lib}"

ls stru/poscar > ls.log
for aa in $(cat ls.log)
do
    rm -fr $aa
    mkdir $aa
    cp INCAR* $aa
    cp vasp.sh $aa
    cp stru/poscar/$aa $aa
    cd $aa
    cp $aa POSCAR
    
    els=`sed -n '6p' POSCAR`
    rm -rf POTCAR
    for el in $els
    do
        if [ $el = 'La' ]; then
            cat $pot_PATH/$el >> POTCAR
        elif [ $el = 'Y' ]; then
            cat $pot_PATH/$el >> POTCAR
        elif [ $el = 'Ce' ]; then
            cat $pot_PATH/$el >> POTCAR
        elif [ $el = 'Th' ]; then
            cat $pot_PATH/$el >> POTCAR
        elif [ $el = 'Be' ]; then
            cat $pot_PATH/$el >> POTCAR
        else
            cat $pot_PATH/H >> POTCAR
        fi
    done
    cd ..
done

