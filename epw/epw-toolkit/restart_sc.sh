#!/bin/bash

sed -e 's/ep_coupling = .true./ep_coupling = .false./' \
    -e 's/elph        = .true./elph        = .false./' \
    -e 's/ephwrite = .true./ephwrite = .false./g'      \
    -e 's/selecqread = .false./selecqread = .true./'   \
    epw_iso_sc.in > new_epw_iso_sc.in

sed -e 's/ep_coupling = .true./ep_coupling = .false./' \
    -e 's/elph        = .true./elph        = .false./' \
    -e 's/ephwrite = .true./ephwrite = .false./g'      \
    -e 's/selecqread = .false./selecqread = .true./'   \
    epw_aniso_sc.in > new_epw_aniso_sc.in

sed -e 's/epw_iso_sc.in/new_epw_iso_sc.in/'            \
    -e 's/epw_iso_sc.out/new_epw_iso_sc.out/'          \
    -e 's/epw_aniso_sc.in/new_epw_aniso_sc.in/'        \
    -e 's/epw_aniso_sc.out/new_epw_aniso_sc.out/'      \
    j5_epw_sc.sh > new_j5_epw_sc.sh