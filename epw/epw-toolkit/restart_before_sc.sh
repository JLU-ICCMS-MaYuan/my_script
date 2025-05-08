#!/bin/bash

prefix=$(grep "prefix" epw_iso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")

echo "Restart from an interrupted q-point while writing ephmatXX files"
echo "Required files: 
1. tmp/${prefix}.epmatwp
2. ${prefix}.ukk
3. crystal.fmt
4. epwdata.fmt
5. vmedata.fmt (or dmedata.fmt)
6. restart.fmt
7. selecq.fmt (selecq.fmt only needed if selecqread = .true. otherwise it will be re-created)
I will link them to the current directory.

This requires to use the same number of pools (npool) as in the original run
"


path=$1

# 判断 path 是否为空
if [ -n "$path" ]; then
    ln -s ${path}/${prefix}.ukk .
    ln -s ${path}/epwdata.fmt .
    ln -s ${path}/vmedata.fmt .
    ln -s ${path}/dmedata.fmt .
fi

sed -i -e 's/ep_coupling = .false./ep_coupling = .true./' \
       -e 's/elph        = .false./elph        = .true./' \
       -e 's/ephwrite = .false./ephwrite = .true./g'      \
       -e 's/selecqread = .true./selecqread = .false./'   \
    epw_iso_sc.in

sed -i -e 's/ep_coupling = .false./ep_coupling = .true./' \
       -e 's/elph        = .false./elph        = .true./' \
       -e 's/ephwrite = .false./ephwrite = .true./g'      \
       -e 's/selecqread = .true./selecqread = .false./'   \
    epw_aniso_sc.in 
