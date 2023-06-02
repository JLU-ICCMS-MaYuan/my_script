#!/bin/sh 

#pot_PATH=/work/home/zhangpy/workplace/HEA/SQS/7/pot
#pot_PATH=/work/home/zhangpy/workplace/HEA/SQS/8-ReW/2/pot

pot_PATH=/public/home/liuhanyu/workplace/mayuan/sqs/prepare-coshare

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
        elif [ $el = 'Be' ]; then
            cat $pot_PATH/$el >> POTCAR
        else
            cat $pot_PATH/H >> POTCAR
        fi
    done
    cd ..
done

