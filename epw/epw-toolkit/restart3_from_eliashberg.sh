#!/bin/bash
prefix=$(grep "prefix" epw_iso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")


echo " Restart by reading ephmatXX files in ${prefix}.ephmat.

This is useful when changing the following input values: nkf1, nkf2, nkf3, nkq1, nkq2, nkq3, fsthick, and degaussw. 

If the message “Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep
to file” has already been output in the previous calculation, the electron-phonon interaction
in the Wannier representation output to ${prefix}.epmatwp can be read in subsequent calculations. 

The ${prefix}.ephmat, restart.fmt, selecq.fmt, and ${prefix}.a2f output in the previous calculation will be useless, you can abandon them or delete them.

Required files: 
1. tmp/${prefix}.ephmat/*   ${prefix}.ephmat directory (which contains egnv, freq, ikmap, ephmatXX files)
2. ${prefix}.dos
3. crystal.fmt
4. selecq.fmt

the path you can input reletive path or absolute path
"


path=$1
abspath=`realpath ${path}`
echo "absolute path : $abspath"

mkdir -p tmp/${prefix}.ephmat
# 判断 path 是否为空
if [ -n "$abspath" ]; then
    ln -sf ${abspath}/tmp/${prefix}.ephmat/* tmp/${prefix}.ephmat/
    ln -sf ${abspath}/tmp/${prefix}.epmatwp  tmp/ 
    ln -sf ${abspath}/${prefix}.ukk .
    ln -sf ${abspath}/crystal.fmt   .
    ln -sf ${abspath}/restart.fmt   .
    ln -sf ${abspath}/selecq.fmt   .
fi

sed -i -e 's/ep_coupling = .false./ep_coupling = .true./'  \
       -e 's/elph        = .false./elph        = .true./'  \
       -e 's/epbwrite    = .true./epbwrite    = .false./'  \
       -e 's/epwwrite    = .true./epwwrite    = .false./'  \
       -e 's/epwread     = .false./epwread     = .true./'  \
       -e 's/wannierize  = .true./wannierize  = .false./'  \
       -e 's/selecqread = .true./selecqread = .false./'    \
       epw_iso_sc.in 


sed -i -e 's/ep_coupling = .false./ep_coupling = .true./'  \
       -e 's/elph        = .false./elph        = .true./'  \
       -e 's/epbwrite    = .true./epbwrite    = .false./'  \
       -e 's/epwwrite    = .true./epwwrite    = .false./'  \
       -e 's/epwread     = .false./epwread     = .true./'  \
       -e 's/wannierize  = .true./wannierize  = .false./'  \
       -e 's/selecqread = .true./selecqread = .false./'    \
       epw_aniso_sc.in 
