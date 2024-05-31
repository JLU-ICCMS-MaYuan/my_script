#!/bin/bash
sed '1,2d' *.dyn0 > qlist.dat
nq=`sed -n 2p *.dyn0`
for Q in `seq 1 $nq`
do
mkdir $Q
cp 1step.sh 2step.sh 3step.sh scf.fit.in scf.in ph.in $Q
Q1=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $1 }')
Q2=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $2 }')
Q3=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $3 }')
cat ./ph.in | sed -e "s/XQ1/$Q1/g" \
| sed -e "s/XQ2/$Q2/g" \
| sed -e "s/XQ3/$Q3/g" \
> $Q/ph.in
done
