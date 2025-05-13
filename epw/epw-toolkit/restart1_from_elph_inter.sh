#!/bin/bash

if [ -f epw_iso_sc.in ]; then
    prefix=$(grep "prefix" epw_iso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_iso_sc.in: $prefix"
elif [ -f epw_aniso_sc.in ]; then
    prefix=$(grep "prefix" epw_aniso_sc.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_aniso_sc.in: $prefix"
elif [ -f epw_elph.out ]; then
    prefix=$(grep "prefix" epw_elph.in | sed -n "s/.*prefix *= *'\(.*\)'.*/\1/p")
    echo "prefix from epw_elph.in: $prefix"
else
    echo "No relevant input or output files found."
fi

echo "Restart by reading ${prefix}.epmatwp and epwdata.fmt.

This is useful when changing the following input values: nkf1, nkf2, nkf3, nkq1, nkq2, nkq3, fsthick, and degaussw. 

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
fi

# sed 替换：无论哪些文件存在，只要有就改
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
