#!/bin/bash

# 使用的代码是MTP结合phonopy
# 需要先安装phonopy以及mlip-2以及MLIP_Phonopy
# MLIP_Phonopy: https://gitlab.com/ivannovikov/mlip_phonopy/-/tree/master
# 编译MLIP_Phonopy时候，最好用intel编译器，gcc编译器编不过去

export PHONOPY_PATH=/work/home/mayuan/miniconda3/bin/
export MLIP_PHONOPY_PATH=/work/home/mayuan/code/mlip_phonopy/

DIM_X=$1
DIM_Y=$2
DIM_Z=$3
"$PHONOPY_PATH"phonopy -d --dim="$DIM_X $DIM_Y $DIM_Z" -c POSCARunitcell

#"$PHONOPY_PATH"phonopy -d --dim="6 6 1" -c POSCARunitcell
mv SPOSCAR POSCAR
num=$(find -type f -name "POSCAR-*" | wc -l)
echo "There are $num supercells"
"$MLIP_PHONOPY_PATH"Main $num
cp POSCARunitcell POSCAR
"$PHONOPY_PATH"phonopy  -p -s band.conf
"$PHONOPY_PATH"phonopy-bandplot  --gnuplot band.yaml > B.txt

grep " frequency" band.yaml  > freq.txt
grep "group_velocity" band.yaml  > gv.txt
