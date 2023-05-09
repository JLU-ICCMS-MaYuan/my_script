#!/usr/bin/env python3

import os
import sys
import shutil

import numpy as np

def change_degauss(old_scf_in_path, new_scf_in_path, degauss):

    with open(old_scf_in_path, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "degauss" in line:
            lines[idx] = f" degauss = {degauss}\n"
        if "pseudo_dir" in line:
            lines[idx] = " pseudo_dir='../pp'\n"

    with open(new_scf_in_path, "w") as f:
        f.writelines(lines)
    
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
#killall -9 pw.x ph.x

mpirun -np 28 /public/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <scf.in> scf.out 
"""

    tangB_slurm = """#!/bin/sh 
#SBATCH  --job-name=test                      
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

mpirun -np 48 /work/home/may/software/qe-7.1/bin/pw.x -npool 4 <scf.in> scf.out 
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

def get_total_energy(scf_out_path):
    """计算并返回每原子的能量 单位的Ry"""
    try:
        Et = float(os.popen(f""" grep "!    total energy " {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        N = int(os.popen(f""" grep "number of atoms/cell" {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        Et_per_atom = Et/N # 注意单位是Ry
        return Et_per_atom
    except:
        print("Try to get total energy from {} failed So let total energy per atom = 1000000000000.0!".format(scf_out_path))   
        Et_per_atom = 1000000000000.0
        return Et_per_atom

if __name__ == "__main__":
    
    # jobsystem = "chaoyin_pbs"
    jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"

    print("Note: --------------------")
    print("    你需要在当前目录下准备好: scf, pp(pp目录中放好赝势)")
    print("    测试的ecutwfc值分别是: 40  50  60  70  80  90")
    print("    测试的ecutrho值分别是: 480 600 720 840 960 10800")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in 0.001 0.01 0.02 0.03 0.04 0.05 0.08 0.1; do cd $i; qsub submit.sh;   cd ..; done")
    print("        for i in 0.001 0.01 0.02 0.03 0.04 0.05 0.08 0.1; do cd $i; sbatch submit.sh; cd ..; done")

    print("Note: --------------------")
    print("    创建测试VASP的ENCUT输入文件目录以及准备vasp的输入文件")
    degausses = [0.1, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001]
    old_scf_in_path = os.path.abspath("scf.in")
    for degauss in degausses:
        test_path = os.path.abspath(str(degauss))
        new_scf_in_path = os.path.join(test_path, "scf.in")
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        change_degauss(old_scf_in_path, new_scf_in_path, degauss)
        write_submit(jobsystem, test_path)
    
    print("    尝试获得每原子的焓值 Et(meV/atom)")
    degauss_Et = []
    for degauss in degausses:
        scf_out_path = os.path.join(str(degauss), "scf.out")
        Et_per_atom = get_total_energy(scf_out_path)
        degauss_Et.append([degauss, Et_per_atom])

    degauss_Et = np.array(degauss_Et)
    Ediff = np.diff(degauss_Et, axis=0)[:,-1]
    Ediff = np.insert(Ediff, 0, [0.0])
    Ediff = Ediff[:, np.newaxis]
    degauss_Et_Ediff = np.hstack((degauss_Et, Ediff))
    if len(degauss_Et) == 8:
        print("{:<12},{:<14},{:<14},{:<14}".format("degauss", "Et(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)"))
        with open("degauss.csv", 'w') as f:
            f.write("{:<12},{:<14},{:<14},{:<14}\n".format("degauss", "Et(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)"))
            for degauss, Et_per_atom, Ediff in degauss_Et_Ediff:
                f.write("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}\n".format(degauss, Et_per_atom, Et_per_atom*13.605693009, Ediff*13.605693009*1000))
                print("{:<12.3f},{:<14.8f},{:<14.8f},{:<14.8f}".format(degauss, Et_per_atom, Et_per_atom*13.605693009, Ediff*13.605693009*1000))
        print("All scf.outs are OK, degauss.csv has been wroten in current position")
    else:
        print("If all OUTCARs are OK, degauss.csv will be wroten in current position")
