#!/bin/bash

# 定义要检查的文件
files=("scffit.in" "scf.in" "s5_PhAssignQ.sh" "split_ph.in")

# 遍历每个文件并检查是否存在
for file in "${files[@]}"; do
    if [[ -f $file ]]; then
        echo "$file exists."
    else
        echo "$file does not exist."
    fi
done

for j in {1..20}; do 
    mkdir $j
    cp inter_dyn_${j} ${j}/Nb4H14.dyn
done

sed '1,2d' Nb4H14.dyn0 > qlist.dat
nq=`sed -n 2p Nb4H14.dyn0`

for Q in `seq 1 $nq`; do
    #mkdir $Q
    cp scffit.in scf.in s5_PhAssignQ.sh $Q
    Q1=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $1 }')
    Q2=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $2 }')
    Q3=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $3 }')
    cat ./split_ph.in | sed -e "s/XQ1/$Q1/g" \
    | sed -e "s/XQ2/$Q2/g" \
    | sed -e "s/XQ3/$Q3/g" \
    > $Q/split_ph.in
done


mkdir 1/tmp
mkdir 1/tmp/_ph0
mkdir 1/tmp/_ph0/Nb4H14.phsave
cp  ../2.fine/1/tmp/_ph0/Nb4H14.Nb4H14.dv*             1/tmp/_ph0
cp  ../2.fine/1/tmp/_ph0/Nb4H14.phsave/patterns.1.xml  1/tmp/_ph0/Nb4H14.phsave

for i in `seq 2 $nq`; do
    cp inter_dyn_$i Nb4H14.dyn
    mkdir $i/tmp
    mkdir $i/tmp/_ph0
    mkdir $i/tmp/_ph0/Nb4H14.phsave
    mkdir $i/tmp/_ph0/Nb4H14.q_${i}
    #cp ../AfterRelax/224/$i/tmp/_ph0/Nb4H14.q_${i}/Nb4H14.Nb4H14.dv* $i/tmp/_ph0/Nb4H14.q_${i}
    cp ../2.fine/$i/tmp/_ph0/Nb4H14.q_${i}/Nb4H14.Nb4H14.dv*    $i/tmp/_ph0/
    cp ../2.fine/$i/tmp/_ph0/Nb4H14.phsave/patterns.${i}.xml    $i/tmp/_ph0/Nb4H14.phsave/patterns.1.xml
    sed -i "4s/$i/1/g" $i/tmp/_ph0/Nb4H14.phsave/patterns.1.xml
done