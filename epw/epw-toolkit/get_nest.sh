#!/bin/bash

band_dat_path=$1

prefix=$(grep "prefix" epw_fermi_nest.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")


echo "# iq Nesting function (q)" > "${prefix}".nesting_fn
grep "Nesting function (q)= " epw_fermi_nest.out | awk '{print $4}' | nl >> "${prefix}".nesting_fn.tmp.dat
paste  "${band_dat_path}" "${prefix}".nesting_fn | awk '{print $1, $4}'  >> "${prefix}".nesting_fn