#!/usr/bin/env python3
import os
import sys
import shutil

import numpy as np

def write_incar(pstress, input_incar, work_path):

    with open(input_incar, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "PSTRESS" in line:
            lines[idx] = f" PSTRESS = {pstress*10}\n"

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

for i in 1 2 3
do 
mpirun -n 28 /public/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1 
cp CONTCAR POSCAR
done 
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


for i in 1 2 3
do 
mpirun -n 48 /work/home/may/software/vasp.6.1.0/bin/vasp_std > vasp.log 2>&1 
cp CONTCAR POSCAR
done 
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

for i in 1 2 3
do 
mpirun -n 64 /public/software/apps/vasp/intelmpi/5.4.4/bin/vasp_std > vasp.log 2>&1
cp CONTCAR POSCAR
done 

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

def get_totalenergy(outcar_path):
    eV2hartree = 0.036749322175655
    hartree2eV = 27.211386245988
    try:
        dE = float(os.popen(f'grep -s "free  energy   TOTEN" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        dE = dE * eV2hartree 
    except:
        dE = 100000000000
    return dE

def get_volume(outcar_path):
    try:
        V = float(os.popen(f'grep -s "volume of cell" {outcar_path}' + " | tail -n 1 | awk '{print $ 5}'").read().strip('\n'))
        V = V / np.power(0.529, 3)
    except:
        V = 100000000000
    return V


if __name__ == "__main__":
    
    jobsystem = "chaoyin_pbs"
    # jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    print("Note: --------------------")
    print("    你需要在当前目录下准备好: POSCAR, POTCAR, INCAR")
    print("    测试的PSTRESS值分别是: 0.0001 50 100 150 200 250 300 350 400")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; qsub submit.sh; cd ..; done")
    print("        for i in p+0.00001 p+50 p+100 p+150 p+200 p+250 p+300 p+350 p+400; do cd $i; sbatch submit.sh; cd ..; done")
    print("\nNote: --------------------")
    print("    创建测试VASP的SIGMA输入文件目录以及准备vasp的输入文件")

    pstress_s = [0.00001, 50, 100, 150, 200, 250, 300, 350, 400]
    potcar_path = os.path.abspath("POTCAR")
    poscar_path = os.path.abspath("POSCAR")
    incar_path  = os.path.abspath("INCAR")

    for pstress in pstress_s:
        if type(pstress) == int :
            test_path = os.path.abspath('p{:+d}'.format(pstress))
        elif type(pstress) == float:
            test_path = os.path.abspath('p{:+.5f}'.format(pstress))
            
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        os.system(f"cp -f {potcar_path} {test_path}")
        os.system(f"cp -f {poscar_path} {test_path}")
        write_incar(pstress, incar_path, test_path)
        write_submit(jobsystem, test_path)

    V_pstress = []
    v_energy = []
    for pstress in pstress_s:
        if type(pstress) == int :
            outcar_path = os.path.join('p{:+d}'.format(pstress), "OUTCAR")
        elif type(pstress) == float:
            outcar_path = os.path.join('p{:+.5f}'.format(pstress), "OUTCAR")
        E = get_totalenergy(outcar_path)
        V = get_volume(outcar_path)
        V_pstress.append([V, pstress])
        v_energy.append([V, E])


    V_pstress = np.array(V_pstress)
    if len(V_pstress) == 9:
        print("{:<14},{:<14}".format("V(A^3)", "pstress(GPa)"))
        with open("V_pstress.csv", 'w') as f:
            f.write("{:<14},{:<14}\n".format("V(A^3)", "pstress(GPa)"))
            for V, pstress in V_pstress:
                f.write("{:<14.8f},{:<14.8f}\n".format(V, pstress))
                print("{:<14.8f},{:<14.8f}".format(V, pstress))
        print("All OUTCARs are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, V_pstress.csv will be wroten in current position")

    v_energy = np.array(v_energy)
    if len(v_energy) == 9:
        print("{:<14},{:<14}".format("V(bohr^3)", "E(hartree)", ))
        with open("eos.in", 'w') as f:
            name             = input("晶体结构名称(建议:crystal)\n");             f.write("{}\n".format(name))
            numatoms         = input("对应原子总数目(建议:胞内原子数)\n");         f.write("{}\n".format(numatoms))
            bm_type          = input("选择的是第几种物态方程(共7种, 建议:3)\n");   f.write("{}\n".format(bm_type))
            insert_points    = input("拟合时每两对E-V之间的插点数目(建议:1000)\n");f.write("{}  {}  {}\n".format(v_energy[0][0], v_energy[-1][0], insert_points))
            num_press_points = input("共有多少个压力点(自己数)\n");               f.write("{}\n".format(num_press_points))
            for V, E in v_energy:
                f.write("{:<14.8f}  {:<14.8f}\n".format(V, E))
                print("{:<14.8f}  {:<14.8f}".format(V, E))
        print("All OUTCARs are OK, V_pstress.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, V_pstress.csv will be wroten in current position")
