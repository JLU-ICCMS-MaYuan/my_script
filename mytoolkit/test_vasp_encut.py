#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np

def write_incar(encut, input_incar, work_path):

    with open(input_incar, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "ENCUT" in line:
            lines[idx] = f" ENCUT = {encut}\n"

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
    print("    测试的ENCUT值分别是: 400, 500, 600, 700, 800, 900, 1000  1100 1200 1300 1400")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的ENCUT输入文件目录以及准备vasp的输入文件")
    encuts = [400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150,]# 1200, 1250, 1300, 1350, 1400, 1450, 1500]
    potcar_path = os.path.abspath("POTCAR")
    poscar_path = os.path.abspath("POSCAR")
    incar_path  = os.path.abspath("INCAR")
    for encut in encuts:
        test_path = os.path.abspath(str(encut))
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        os.system(f"cp -f {potcar_path} {test_path}")
        os.system(f"cp -f {poscar_path} {test_path}")
        write_incar(encut, incar_path, test_path)
        write_submit(jobsystem, test_path)

    print("    尝试获得每原子的焓值 dE(meV/atom)")
    encut_dE_dF_dV = []
    for encut in encuts:
        outcar_path = os.path.join(str(encut), "OUTCAR")
        dE_per_atom, force, virial_per_atom = get_info_per_atom(outcar_path)
        encut_dE_dF_dV.append([encut, dE_per_atom, force, virial_per_atom])


    # 方式一：通过相邻两个encut的能量的差判断收敛性
    # Ediff = np.diff(encut_dE, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # encut_dE_Hdiff = np.hstack((encut_dE, Ediff))

    # 方式二：所有encut的能量减去最高encut的能量判断收敛性
    with open("encut_dE.csv", 'w') as f:
        f.write("{:<12},{:<14},{:<14},{:<14}\n".format("encut", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        print("{:<12},{:<14},{:<14},{:<14}".format("encut", "Ediff(meV/atom)", "Fdiff(meV/A)", "Vdiff(meV/atom)"))
        for idx, eefv in enumerate(encut_dE_dF_dV):
            kmeshes = eefv[0] 
            delta_E = eefv[1] - encut_dE_dF_dV[-1][1]
            delta_F = eefv[2] - encut_dE_dF_dV[-1][2]
            delta_V = eefv[3] - encut_dE_dF_dV[-1][3]
            f.write("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}\n".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000,np.nanmax(delta_V)*1000))
            print("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}".format(kmeshes, delta_E*1000, np.nanmax(delta_F)*1000, np.nanmax(delta_V)*1000))
        else:
            print("If all OUTCARs are OK, encut_dE.csv will be wroten in current position")

