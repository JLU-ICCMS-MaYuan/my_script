#!/bin/bash

for aa in `cat ls.log` ; do 
    cd $aa
    sbatch vasp.sh
    cd ..
done
