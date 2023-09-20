#!/bin/bash

for aa in `cat ls.log` ; do 
    cd $aa
    qsub vasp.sh
    cd ..
done
