#!/bin/bash

echo "This scripts had better be executed in 1.dp-data/2.testset"

ln -s ../1.trainset/* . 
iterpaths=`find "../../2.getdp-mod" -maxdepth 1 -type d -name "iter*"`
for iterpath in $iterpaths; do
       # echo "$iterpath"
       datapaths=`find "${iterpath}/02.fp" -type d -name "data*" 2>/dev/null`
       for datapath in $datapaths; do
               # echo "$datapath"
               dataname=$(basename "$datapath")
               itername=$(basename "$iterpath")
               newname="${itername}_${dataname}"
               echo "$newname"
               ln -s "$datapath" "$newname"
       done
done
