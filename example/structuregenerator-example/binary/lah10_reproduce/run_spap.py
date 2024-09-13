#!/usr/bin/python3

from spap import run_spap

run_spap(
    total_struc=150,
    # l_cif=True,
    # r_cut_off=7.5,
    symprec=0.02,
    l_poscar=True,
    threshold=0.045,
    # l_comp=False,
    # ilat=1,
    work_dir='./',
    i_mode=2,
)

