#!/bin/bash
killall -9 vasp_std
#mpirun -n 64 /public/software/apps/vasp/intelmpi/6.3.2/bin/vasp_std > vasp.log 2>&1
#mpirun -n 48 /work/software/vasp.6.1.0/vasp_std > vasp.log 2>&1
mpirun -np 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
