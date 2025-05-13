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
        awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $3*1000}' ${iso_sc_dir}/${prefix}.imag_iso_* > ${prefix}.IMAG_ISO_MUC_${mu}_GAP0
    else
        echo "Warning: Directory $sc_dir not found. Trying current directory..."
        awk 'FNR==2 {n=split(FILENAME,a,"_"); temp=a[length(a)]; print FILENAME, temp, $3*1000}' ${prefix}.imag_iso_* > ${prefix}.IMAG_ISO_MUC_${mu}_GAP0
        if [ -s "${prefix}.IMAG_ISO_MUC_${mu}_GAP0" ]; then
            echo "Data found in ${iso_sc_dir}    Stopping script."
            cat ${prefix}.IMAG_ISO_MUC_${mu}_GAP0
            exit 0
        fi

    echo "aniso_muc_${mu}"
    aniso_sc_dir="aniso_muc_${mu}"
    if [ -d "$aniso_sc_dir" ]; then
        cat ${aniso_sc_dir}/${prefix}.imag_aniso_gap0_* > ${prefix}.IMAG_ANISO_MUC_${mu}_GAP0
    else
        echo "Warning: Directory $sc_dir not found. Trying current directory..."
        cat ${prefix}.imag_aniso_gap0_* > ${prefix}.IMAG_ANISO_MUC_${mu}_GAP0
        if [ -s "${prefix}.IMAG_ANISO_MUC_${mu}_GAP0" ]; then
            echo "Data found in ${aniso_sc_dir}    Stopping script."
            cat ${prefix}.IMAG_ANISO_MUC_${mu}_GAP0
            exit 0
        fi
done