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

def get_energy_per_atom(outcar_path):
    """计算并返回每原子的焓值, 单位是eV"""
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
    print("    测试的ENCUT值分别是: 400, 500, 600, 700, 800, 900, 1000  1100 1200 1300 1400")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 400 500 600 700 800 900 1000 1100 1200 1300 1400; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 400 500 600 700 800 900 1000 1100 1200 1300 1400; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的ENCUT输入文件目录以及准备vasp的输入文件")
    encuts = [400, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500]
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
    encut_dE = []
    for encut in encuts:
        outcar_path = os.path.join(str(encut), "OUTCAR")
        dE_per_atom = get_energy_per_atom(outcar_path)
        encut_dE.append([encut, dE_per_atom])

    encut_dE = np.array(encut_dE)

    # 方式一：通过相邻两个encut的能量的差判断收敛性
    # Ediff = np.diff(encut_dE, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # encut_dE_Hdiff = np.hstack((encut_dE, Ediff))

    # 方式二：所有encut的能量减去最高encut的能量判断收敛性
    Ediff = encut_dE[:, 1] - encut_dE[-1, -1]
    Ediff = Ediff[:, np.newaxis]
    encut_dE_Hdiff = np.hstack((encut_dE, Ediff))
    
    if len(encut_dE_Hdiff) == 21:
        print("{:<12},{:<14},{:<14}".format("ENCUT", "dE(eV/atom)", "diff(meV/atom)"))
        with open("encut_dE.csv", 'w') as f:
            f.write("{:<12},{:<14},{:<14}\n".format("ENCUT", "dE(eV/atom)", "diff(meV/atom)"))
            for encut, dE_per_atom, Hdiff in encut_dE_Hdiff:
                f.write("{:<12.3f},{:<14.8f},{:<14.8f}\n".format(encut, dE_per_atom, Hdiff*1000))
                print("{:<12.3f},{:<14.8f},{:<14.8f}".format(encut, dE_per_atom, Hdiff*1000))
        print("All OUTCARs are OK, encut_dE.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, encut_dE.csv will be wroten in current position")

