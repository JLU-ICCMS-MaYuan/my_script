#!/bin/sh
# set the needed directories and files
PREFIX=`grep prefix scf.in | cut -d "'" -f 2`
nq=`sed -n '2p' *dyn0`
#rm -rf elph_dir
mkdir elph_dir
for n in $(seq 1 $nq)
do
cd $n/elph_dir
cp ../$PREFIX.dyn ../../$PREFIX.dyn$n
cp elph.inp_lambda.1 ../../elph_dir/elph.inp_lambda.$n
cp a2Fq2r.51.1 ../../elph_dir/a2Fq2r.51.$n
cp a2Fq2r.52.1 ../../elph_dir/a2Fq2r.52.$n
cp a2Fq2r.53.1 ../../elph_dir/a2Fq2r.53.$n
cp a2Fq2r.54.1 ../../elph_dir/a2Fq2r.54.$n
cp a2Fq2r.55.1 ../../elph_dir/a2Fq2r.55.$n
cp a2Fq2r.56.1 ../../elph_dir/a2Fq2r.56.$n
cp a2Fq2r.57.1 ../../elph_dir/a2Fq2r.57.$n
cp a2Fq2r.58.1 ../../elph_dir/a2Fq2r.58.$n
cp a2Fq2r.59.1 ../../elph_dir/a2Fq2r.59.$n
cp a2Fq2r.60.1 ../../elph_dir/a2Fq2r.60.$n
echo ../../elph_dir/elph.inp_lambda.$n
cd ../..
done
