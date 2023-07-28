#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np

def write_incar(kspacing, input_incar, work_path):

    with open(input_incar, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "KSPACING" in line:
            lines[idx] = f" KSPACING = {kspacing}\n"

    incar_path = os.path.join(work_path, "INCAR")
    with open(incar_path, "w") as incarfile:
        incarfile.writelines(lines)

def write_submit(jobtype:str, work_path):

    chaoyin_pbs = """#!/bin/sh
#PBS -N    mayqe                                    
#PBS -q    liuhy         
#PBS -l    nodes=1:ppn=28               
#PBS -j    oe                                      
#PBS -V  
source /public/home/mayuan/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited        
cd $PBS_O_WORKDIR    

#killall -9 vasp_std

mpirun -n 28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1  
"""

    tangB_slurm = """#!/bin/sh
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=lhy          
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         

source /work/home/may/intel/oneapi/setvars.sh --force
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
ulimit -s unlimited

mpirun -n 48 /work/home/may/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1 
"""

    coshare_slurm = """#!/bin/sh
#SBATCH  --job-name=vaspopt-1
#SBATCH  --output=log.out
#SBATCH  --error=log.err
#SBATCH  --partition=normal                       
#SBATCH  --nodes=1                            
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                    
#SBATCH  --exclude=node37,node34,node38,node48,node10,node15,node5,node20

source /public/env/mpi_intelmpi-2021.3.0.sh
source /public/env/compiler_intel-compiler-2021.3.0.sh

ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
export I_MPI_FABRICS=shm
export MKL_DEBUG_CPU_TYPE=5

mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
"""

    if jobtype == "chaoyin_pbs":
        jobsystem = chaoyin_pbs
    elif jobtype == "tangB_slurm":
        jobsystem = tangB_slurm
    elif jobtype == "coshare_slurm":
        jobsystem = coshare_slurm
    else:
        print("其它机器的任务系统, 你没有准备这个机器的提作业脚本. 程序退出. ")
        sys.exit(1)

    submit_path = os.path.join(work_path, "submit.sh")
    with open(submit_path, "w") as jobfile:
        jobfile.write(jobsystem)

def get_energy_per_atom(outcar_path):
    try:
        dE = float(os.popen(f'grep -s "free  energy   TOTEN" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        begin_id = os.popen(f'grep -n "position of ions in cartesian coordinates" {outcar_path}').read().split(":")[0]
        N = 0; row_id=int(begin_id)
        while True:
            row_id = row_id+1
            content  = os.popen(f"sed -n '{row_id}p' {outcar_path}").read().strip().split()
            if len(content) == 3:
                N += 1
            else:
                break
        dE_per_atom = dE/N
        return dE_per_atom
    except:
        print("Try to get DeltaH from {} failed So let dE_per_atoms = 1000000000000.0!".format(outcar_path))   
        dE_per_atom = 1000000000000.0
        return dE_per_atom



if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, POTCAR, INCAR")
    print("    测试的KSPACING值分别是: 0.3, 0.2, 0.18, 0.15, 0.12")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 0.4 0.3 0.25 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 0.4 0.3 0.25 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的KSPACING输入文件目录以及准备vasp的输入文件")
    kspacings = [0.4, 0.3, 0.25, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12]
    potcar_path = os.path.abspath("POTCAR")
    poscar_path = os.path.abspath("POSCAR")
    incar_path  = os.path.abspath("INCAR")

    for kspacing in kspacings:
        test_path = os.path.abspath(str(kspacing))
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        
        os.system(f"cp -f {potcar_path} {test_path}")
        os.system(f"cp -f {poscar_path} {test_path}")
        write_incar(kspacing, incar_path, test_path)
        write_submit(jobsystem, test_path)

    print("    尝试获得每原子的焓值 dE(meV/atom)")
    kspacing_dE = []
    for kspacing in kspacings:
        outcar_path = os.path.join(str(kspacing), "OUTCAR")
        dE_per_atom = get_energy_per_atom(outcar_path)
        kspacing_dE.append([kspacing, dE_per_atom])

    kspacing_dE = np.array(kspacing_dE)
    Ediff = np.diff(kspacing_dE, axis=0)[:,-1]
    Ediff = np.insert(Ediff, 0, [0.0])
    Ediff = Ediff[:, np.newaxis]
    kspacing_dE_Hdiff = np.hstack((kspacing_dE, Ediff))
    if len(kspacing_dE_Hdiff) == 12:
        print("{:<12},{:<14},{:<14}".format("kspacing", "dE(eV/atom)", "diff(meV/atom)"))
        with open("kspacing_dE.csv", 'w') as f:
            f.write("{:<12},{:<14},{:<14}\n".format("kspacing", "dE(eV/atom)", "diff(meV/atom)"))
            for kspacing, dE_per_atom, Hdiff in kspacing_dE_Hdiff:
                f.write("{:<12.3f},{:<14.8f},{:<14.8f}\n".format(kspacing, dE_per_atom, Hdiff*1000))
                print("{:<12.3f},{:<14.8f},{:<14.8f}".format(kspacing, dE_per_atom, Hdiff*1000))
        print("All OUTCARs are OK, kspacing_dE.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, kspacing_dE.csv will be wroten in current position")

