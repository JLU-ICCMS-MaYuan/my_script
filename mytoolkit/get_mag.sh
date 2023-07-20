#!/bin/bash

echo "You need to tell me: How many lines do you want to print out? For example: get_meg.sh 20"
echo ""

lineid=$(grep -n "magnetization (x)" OUTCAR | tail -n 1 | awk -F':' '{print $1}')

awk "NR >= $lineid && NR < $lineid + $1" OUTCAR 
