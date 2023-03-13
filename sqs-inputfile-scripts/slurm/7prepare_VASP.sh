#!/bin/sh 

#pot_PATH=/work/home/zhangpy/workplace/HEA/SQS/7/pot
#pot_PATH=/work/home/zhangpy/workplace/HEA/SQS/8-ReW/2/pot

pot_PATH=/work/home/may/POT/vasp_pot1/potpaw_PBE54
ls stru/poscar > ls.log
for aa in $(cat ls.log)
do
    rm -fr $aa
    mkdir $aa
    cp INCAR $aa
    cp vasp.sh $aa
    cp stru/poscar/$aa $aa
    cd $aa
    cp $aa POSCAR
    
    els=`sed -n '6p' POSCAR`
    rm -rf POTCAR
    for el in $els
    do
        if [ $el = 'La' ]; then
            cat $pot_PATH/$el/POTCAR >> POTCAR
        elif echo "Ba_sv" | grep -q "$el"; then
            cat $pot_PATH/Ba_sv/POTCAR >> POTCAR
        else
            cat $pot_PATH/H/POTCAR >> POTCAR
        fi
    done
    cd ..
done

