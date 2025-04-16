#!/bin/bash

#accroding to the symmetry of system, creat the disorder structure.
# -2 means the distance of two cluster formed by two atoms. Usually, -2 is the value between nearest and second-nearest neighbons.
# -3 means the distance of two cluster formed by three atoms.

corrdump -l=rndstr.in -ro -noe -nop -clus -2=4.0 -3=3.6 -4=3.5; getclus

