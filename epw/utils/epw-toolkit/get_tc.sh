#!/bin/bash

if [ -n "$1" ]; then
    prefix=$1
    echo "prefix from command-line argument: $prefix"
else
    # 自动从输入文件中提取 prefix
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
        echo "No relevant input or output files found. Please provide a prefix."
        exit 1
    fi
fi


for mu in 0.1 0.13 0.16; do

    echo "iso_muc_${mu}"
    iso_sc_dir="iso_muc_${mu}"
    if [ -d "$iso_sc_dir" ]; then
        awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $3*1000}' ${iso_sc_dir}/${prefix}.imag_iso_* > ${prefix}.IMAG_ISO_MUC_${mu}_GAP0.dat
        awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $4*1000}' ${iso_sc_dir}/${prefix}.pade_iso_* > ${prefix}.PADE_ISO_MUC_${mu}_GAP0.dat
        awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $4*1000}' ${iso_sc_dir}/${prefix}.acon_iso_* > ${prefix}.ACON_ISO_MUC_${mu}_GAP0.dat
    else
        echo "Warning: Directory $sc_dir not found. Trying current directory..."
        if [ -f "epw_iso_sc.in" ]; then
            new_mu=$(grep "muc" epw_iso_sc.in | sed -n "s/.*muc *= *\([0-9.eE+-]*\).*/\1/p")
            awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $3*1000}' ${prefix}.imag_iso_* > ${prefix}.IMAG_ISO_MUC_${new_mu}_GAP0.dat
            awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $4*1000}' ${prefix}.pade_iso_* > ${prefix}.PADE_ISO_MUC_${new_mu}_GAP0.dat
            awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $4*1000}' ${prefix}.acon_iso_* > ${prefix}.ACON_ISO_MUC_${new_mu}_GAP0.dat
         fi
     fi

    echo "aniso_muc_${mu}"
    aniso_sc_dir="aniso_muc_${mu}"
    if [ -d "$aniso_sc_dir" ]; then
        find "${aniso_sc_dir}" -maxdepth 1 -type f -name "${prefix}.imag_aniso_gap0_*" ! -name '*.cube' ! -name '*.frmsf' -exec cat {} + > ${prefix}.IMAG_ANISO_MUC_${mu}_GAP0.dat
        find "${aniso_sc_dir}" -maxdepth 1 -type f -name "${prefix}.pade_aniso_gap0_*"  -exec cat {} + > ${prefix}.PADE_ANISO_MUC_${mu}_GAP0.dat
        find "${aniso_sc_dir}" -maxdepth 1 -type f -name "${prefix}.acon_aniso_gap0_*"  -exec cat {} + > ${prefix}.ACON_ANISO_MUC_${mu}_GAP0.dat
    else
        if [ -f "epw_aniso_sc.in" ]; then
            echo "Warning: Directory $sc_dir not found. Trying current directory..."
            new_mu=$(grep "muc" epw_aniso_sc.in | sed -n "s/.*muc *= *\([0-9.eE+-]*\).*/\1/p")
            find . -maxdepth 1 -type f -name "${prefix}.imag_aniso_gap0_*" ! -name '*.cube' ! -name '*.frmsf' -exec cat {} + > ${prefix}.IMAG_ANISO_MUC_${new_mu}_GAP0.dat
            find . -maxdepth 1 -type f -name "${prefix}.pade_aniso_gap0_*"  -exec cat {} + > ${prefix}.PADE_ANISO_MUC_${new_mu}_GAP0.dat
            find . -maxdepth 1 -type f -name "${prefix}.acon_aniso_gap0_*"  -exec cat {} + > ${prefix}.ACON_ANISO_MUC_${new_mu}_GAP0.dat
        fi
     fi
done


for mu in 0.1 0.13 0.16; do
    for type in IMAG PADE ACON; do
        iso_file="${prefix}.${type}_ISO_MUC_${mu}_GAP0.dat"
        if [ -s "$iso_file" ]; then
            echo "$iso_file is OK!"
            sed -i "1iFilename Temperature(K) ISO_${type}_MUC_${mu}_delta_nk[meV]" "$iso_file"
        else
            echo "$iso_file is empty or does not exist. Deleting..."
            rm -f "$iso_file"
        fi

        aniso_file="${prefix}.${type}_ANISO_MUC_${mu}_GAP0.dat"
        if [ -s "$aniso_file" ]; then
            sed -i "1cTemperature(K) ANISO_${type}_MUC_${mu}_delta_nk[meV]" "$aniso_file"
            echo "$aniso_file is OK!"
        else
            echo "$aniso_file is empty or does not exist. Deleting..."
            rm -f "$aniso_file"
        fi
    done
done
