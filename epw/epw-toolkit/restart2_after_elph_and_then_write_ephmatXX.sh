#!/bin/bash

if [ -f epw_iso_sc.in ]; then
    prefix=$(grep "prefix" epw_iso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_iso_sc.in: $prefix"
elif [ -f epw_aniso_sc.in ]; then
    prefix=$(grep "prefix" epw_aniso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_aniso_sc.in: $prefix"
elif [ -f epw_elph.in ]; then
    prefix=$(grep "prefix" epw_elph.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_elph.in: $prefix"
else
    echo "No relevant input or output files found."
fi

echo "Restart from an interrupted q-point while writing ephmatXX files

Required files: 
1. tmp/${prefix}.epmatwp
2. ${prefix}.ukk
3. crystal.fmt
4. epwdata.fmt
5. vmedata.fmt (or dmedata.fmt)
6. restart.fmt
7. selecq.fmt (selecq.fmt only needed if selecqread = .true. otherwise it will be re-created)
I will link them to the current directory.

This requires to use the same number of pools (npool) as in the original run

the path you can input reletive path or absolute path
"


path=$1
abspath=`realpath ${path}`
echo "absolute path : $abspath"

# 判断 path 是否为空
mkdir -p tmp
if [ -n "$abspath" ]; then
    ln -sf ${abspath}/tmp/${prefix}.epmatwp tmp/
    ln -sf ${abspath}/${prefix}.ukk .
    ln -sf ${abspath}/crystal.fmt   .
    ln -sf ${abspath}/epwdata.fmt   .
    ln -sf ${abspath}/vmedata.fmt   .
    ln -sf ${abspath}/dmedata.fmt   .
    ln -sf ${abspath}/restart.fmt   .
    ln -sf ${abspath}/selecq.fmt    .
fi

for infile in epw_iso_sc.in epw_aniso_sc.in epw_elph.in; do
    if [ -f "$infile" ]; then
    sed -i -e 's/ep_coupling = .false./ep_coupling = .true./'  \
           -e 's/elph        = .false./elph        = .true./'  \
           -e 's/epbwrite    = .true./epbwrite    = .false./'  \
           -e 's/epwwrite    = .true./epwwrite    = .false./'  \
           -e 's/epwread     = .false./epwread     = .true./'  \
           -e 's/wannierize  = .true./wannierize  = .false./'  \
           -e 's/selecqread = .true./selecqread = .false./'    \
               "$infile"
    fi
done
