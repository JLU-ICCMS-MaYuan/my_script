#!/bin/bash


echo "./get_dyn.sh <sparse_qmesh> <sparse_dyns> <fine_qmesh> <fine_dyns> <prefix>"

rm *.dyn*  *_dyn_*

sparse_qmesh=$1
sparse_dyns=$2
fine_qmesh=$3
fine_dyns=$4
prefix=$5

for i in $(seq 1 $sparse_dyns); do
	cp ../1.sparse/${prefix}.dyn$i ${sparse_qmesh}.dyn$i
done

for i in $(seq 1 $fine_dyns); do
	cp ../2.fine/${prefix}.dyn$i ${fine_qmesh}.dyn$i
done

cp ../V3_Hessian.dyn* .
