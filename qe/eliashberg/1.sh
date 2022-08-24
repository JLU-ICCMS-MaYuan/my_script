#!/bin/sh
   for a in eliashberg.f90  init0.f90     modmain.f90  modphonon.f90  pade.f90\
 eliashberg.x    init1.f90     modmain.mod  modphonon.mod  readalpha2f.f90\
 fderiv.f90      install.sh    modmpi.f90   mpi.mod        spline.f90\
 flushifc.f90    mcmillan.f90  modmpi.mod   mpi_stub.f90   timesec.f90
   do 
   cd /work/home/liull/sunyao/EPC/src
   chmod 777 $a
   done
