#!/bin/bash
prefix=$(grep "prefix" epw_iso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")

echo "Restart by reading ${prefix}.epmatwp and epwdata.fmt."

echo "This is useful when changing the following input values: nkf1, nkf2, nkf3, nkq1, nkq2, nkq3, fsthick, and degaussw. 

If the message “Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep
to file” has already been output in the previous calculation, the electron-phonon interaction
in the Wannier representation output to ${prefix}.epmatwp can be read in subsequent calculations. 

The ${prefix}.ephmat, restart.fmt, selecq.fmt, and ${prefix}.a2f output in the previous calculation will be useless, you can abandon them or delete them.

Required files: 
1. tmp/${prefix}.epmatwp
2. ${prefix}.ukk
3. crystal.fmt
4. epwdata.fmt
5. vmedata.fmt (or dmedata.fmt)
"


path=$1

# 判断 path 是否为空
if [ -n "$path" ]; then
    ln -s ${path}/tmp/${prefix}.epmatwp . 
    ln -s ${path}/tmp/${prefix}.prefix.ukk .
    ln -s ${path}/crystal.fmt .
    ln -s ${path}/restart.fmt .
    ln -s ${path}/selectq.fmt .
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
