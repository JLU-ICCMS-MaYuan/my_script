#!/bin/bash
#mpirun -n 64 /public/software/apps/vasp/intelmpi/6.3.2/bin/vasp_std > vasp.log 2>&1
mpirun -np 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
