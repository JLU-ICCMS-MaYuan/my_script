#bin/sh -v
#source /data/home/mym/soft/intel_2019/compilers_and_libraries_2019.4.243/linux/bin/compilervars.sh intel64
#source /data/home/mym/soft/intel_2019/compilers_and_libraries_2019.4.243/linux/mkl/bin/mklvars.sh intel64
#source /data/home/mym/soft/intel_2019/impi/2019.4.243/intel64/bin/mpivars.sh

source ~/intel/oneapi/setvars.sh --force

ifort -o eliashberg.x spline.f90 timesec.f90 fderiv.f90 flushifc.f90 init0.f90 mpi_stub.f90 mcmillan.f90 modmain.f90 modphonon.f90 pade.f90 modmpi.f90 readalpha2f.f90 eliashberg.f90 -O3 -ip -unroll -no-prec-div -parallel
