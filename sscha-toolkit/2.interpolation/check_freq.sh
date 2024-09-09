#!/bin/bash


for i in {1..6}; do echo $i; grep freq inter_dyn_$i | head ; done
