#!/bin/bash

inter_dyns=$1
for i in $(seq 1 $inter_dyns); do echo $i; grep freq inter_dyn_$i | head ; done
