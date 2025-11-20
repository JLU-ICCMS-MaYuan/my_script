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

abspath=`pwd`
for mu in 0.1 0.13 0.16; do
    for sc in iso aniso; do
        mkdir -p ${sc}_${mu}/tmp/${prefix}.ephmat
        echo "${sc}_${mu}"
        cp epw_${sc}_sc.in ${sc}_${mu}
        cp j5_epw_sc.sh    ${sc}_${mu}
        sed -i -e 's/ep_coupling = .true./ep_coupling = .false./'  \
           -e 's/elph        = .true./elph        = .false./'  \
           -e 's/epbwrite    = .true./epbwrite    = .false./'  \
           -e 's/epwwrite    = .true./epwwrite    = .false./'  \
           -e 's/epwread     = .false./epwread     = .true./'  \
           -e 's/ephwrite = .true./ephwrite = .false./'        \
           -e 's/wannierize  = .true./wannierize  = .false./'  \
           -e 's/selecqread = .true./selecqread = .false./'    \
           -e "s/.*mu.*/mu=${mu}/"                             \
        epw_${sc}_sc.in
        ln -sf ${abspath}/tmp/${prefix}.ephmat/* ${sc}_${mu}/tmp/${prefix}.ephmat/
        ln -sf ${abspath}/tmp/${prefix}.epmatwp  ${sc}_${mu}/tmp/ 
        ln -sf ${abspath}/${prefix}.ukk          ${sc}_${mu}
        ln -sf ${abspath}/crystal.fmt            ${sc}_${mu}
        ln -sf ${abspath}/restart.fmt            ${sc}_${mu}
        ln -sf ${abspath}/selecq.fmt             ${sc}_${mu}
     done
done
