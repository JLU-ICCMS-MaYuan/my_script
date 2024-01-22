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

def get_info_per_atom(outcar_path):

    try:
        natoms = int(os.popen("grep 'NIONS' %s |awk '{print $12}'"%(outcar_path)).read())
    except:
        natoms = None

    try:
        energy = float(os.popen(f'grep -s "energy  without entropy" {outcar_path}' + " | awk '{print $4}'").read().strip('\n'))
    except:
        energy = None

    try:
        force = np.array(list(map(float, os.popen("grep -A%d 'TOTAL-FORCE' %s |tail -n %d|awk '{print $4,$5,$6}'"%(natoms+1,outcar_path,natoms)).read().split())))
    except:
        force = None

    try:
        virial = np.array(list(map(float, os.popen("grep -A20 '\-STRESS' %s |grep Total|awk '{print $2,$3,$4,$5,$6,$7}'"%(outcar_path)).read().split())))
    except:
        virial = None
  
    if (natoms is not None) and (energy is not None) and (force is not None) and (virial is not None):
        dE_per_atom = np.array(energy)/natoms
        virial_per_atom = np.array(virial)/natoms
    else:
        dE_per_atom = 10000000000000000
        force = 10000000000000000 
        virial_per_atom = 10000000000000000
    return dE_per_atom, force, virial_per_atom



if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, POTCAR, INCAR")
    print("    测试的KSPACING值分别是: 0.3, 0.2, 0.18, 0.15, 0.12")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 0.4 0.3 0.25 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.1 0.09 0.08; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 0.4 0.3 0.25 0.2 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.1 0.09 0.08; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的KSPACING输入文件目录以及准备vasp的输入文件")
    kspacings = [0.4, 0.3, 0.25, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08]
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
    kspacing_dE_dF_dV = []
    for kspacing in kspacings:
        outcar_path = os.path.join(str(kspacing), "OUTCAR")
        dE_per_atom, force, virial_per_atom = get_info_per_atom(outcar_path)
        kspacing_dE_dF_dV.append([kspacing, dE_per_atom, force, virial_per_atom])


    # 方式一：通过相邻两个kspacing的能量的差判断收敛性
    # Ediff = np.diff(kspacing_dE, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # kspacing_dE_Hdiff = np.hstack((kspacing_dE, Ediff))

    # 方式二：所有kpoints的能量减去最密kspacing的能量判断收敛性

    with open("kspacing_dE.csv", 'w') as f:
        f.write("{:<12},{:<14},{:<14},{:<14}\n".format("kspacing", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        print("{:<12},{:<14},{:<14},{:<14}".format("kspacing", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        for idx, kefv in enumerate(kspacing_dE_dF_dV):
            kmeshes = kefv[0] 
            delta_E = kefv[1] - kspacing_dE_dF_dV[-1][1]
            delta_F = kefv[2] - kspacing_dE_dF_dV[-1][2]
            delta_V = kefv[3] - kspacing_dE_dF_dV[-1][3]
            f.write("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}\n".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000,np.nanmax(delta_V)*1000))
            print("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000, np.nanmax(delta_V)*1000))
        else:
            print("If all OUTCARs are OK, kspacing_dE.csv will be wroten in current position")
