from numpy import *
import numpy as np

import os

n_sub = 100
n_pw=10
#all_scf_files = [os.path.join("run_calculation", f) for f in os.listdir("run_calculation") if f.startswith("espresso_run_")]
header="""#!/bin/sh 
#------ slurm option --------#
#SBATCH --partition=mpc
#SBATCH --account=hp240139
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00

ulimit -s unlimited

#------- Program execution -------#
module load intelmpi/impi_23.2.0
module load intel/23.02.1
"""
for i in range(1,n_sub+1):
    filename= "sub_{}.sh".format(i)
    run_lines=[]
    j_start=(i-1)*n_pw+1
    j_end=i*n_pw
    for j in range(j_start,j_end+1):
        run_line="timeout 1800s srun /home/h240012/soft/qe-7.0/bin/pw.x -nk 2 -in  espresso_run_{}.pwi > espresso_run_{}.pwo 2>&1".format(j,j)
        run_lines.append(run_line)
    with open(filename, "w") as f:
        f.write(header)
        for run_line in run_lines:
            print(run_line,file=f)
    os.system(rf"mv {filename} run_calculation")


