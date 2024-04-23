#!/usr/bin/env python3

import os
import sys

import numpy as np

def create_kmesh(kresolution, old_scf_in_path):
    with open(old_scf_in_path, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "CELL_PARAMETERS" in line:
            begin_id = idx

    ax,ay,az = [ float(x) for x in lines[begin_id+1].strip('\n').split() ]
    bx,by,bz = [ float(x) for x in lines[begin_id+2].strip('\n').split() ]
    cx,cy,cz = [ float(x) for x in lines[begin_id+3].strip('\n').split() ]

    a1 = np.array([ax,ay,az])
    a2 = np.array([bx,by,bz])
    a3 = np.array([cx,cy,cz])

    b1 = 2*np.pi*np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
    b2 = 2*np.pi*np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
    b3 = 2*np.pi*np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))

    b1mol=np.linalg.norm(b1)
    b2mol=np.linalg.norm(b2)
    b3mol=np.linalg.norm(b3)

    n_k1 = np.ceil(b1mol/kresolution)
    n_k2 = np.ceil(b2mol/kresolution)
    n_k3 = np.ceil(b3mol/kresolution)
    print("    kresolution={:<4.3f}, 相应的nk1, nk2, nk3 = {:<4} {:<4} {:<4}".format(kresolution, n_k1, n_k2, n_k3))

    return n_k1, n_k2, n_k3

def change_kpoints(old_scf_in_path, new_scf_in_path, kresolution):

    with open(old_scf_in_path, "r") as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if "pseudo_dir" in line:
            lines[idx] = " pseudo_dir='../pp'\n"

    n_k1, n_k2, n_k3 = create_kmesh(kresolution, old_scf_in_path)
    lines[-1] = str(n_k1) +" "+ str(n_k2) +" "+ str(n_k3) + " 0 0 0" '\n'

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

mpirun -np 64 /public/software/apps/quantum-espresso/intelmpi/6.7/bin/pw.x -npool 8 <scf.in> scf.out
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

    riken_slurm = """#!/bin/bash
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

srun /home/h240012/soft/qe-7.0/bin/pw.x -nk 4 -in scf.in > scf.out 2>&1
    """

    if jobtype == "chaoyin_pbs":
        jobsystem = chaoyin_pbs
    elif jobtype == "tangB_slurm":
        jobsystem = tangB_slurm
    elif jobtype == "coshare_slurm":
        jobsystem = coshare_slurm
    elif jobtype == "riken_slurm":
        jobsystem = riken_slurm
    else:
        print("其它机器的任务系统, 你没有准备这个机器的提作业脚本. 程序退出. ")
        sys.exit(1)

    submit_path = os.path.join(work_path, "submit.sh")
    with open(submit_path, "w") as jobfile:
        jobfile.write(jobsystem)

def get_total_energy(scf_out_path):
    """计算并返回每原子的能量 单位的Ry"""
    try:
        Et = float(os.popen(f""" grep -s "!    total energy " {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        N = int(os.popen(f""" grep -s "number of atoms/cell" {scf_out_path}""" + """  | awk '{print $5}' """).read().strip('\n'))
        Et_per_atom = Et/N # 注意单位是Ry
        return Et_per_atom
    except:
        print("Try to get total energy from {} failed So let total energy per atom = 1000000000000.0!".format(scf_out_path))
        Et_per_atom = 1000000000000.0
        return Et_per_atom

def checktime(scf_out_path):

    try:
        timetot = os.popen(f""" grep PWSCF {scf_out_path} | grep CPU | tail -n 1 |""" + """awk -F':' '{print $2}'""").read().strip('\n')
    except:
        timetot = None

    return timetot

if __name__ == "__main__":

    # jobsystem = "chaoyin_pbs"
    # jobsystem = "tangB_slurm"
    # jobsystem = "coshare_slurm"
    jobsystem = "riken_slurm"

    print("Note: --------------------")
    print("    你需要在当前目录下准备好: scf, pp(pp目录中放好赝势)")

    print("Note: --------------------")
    print("    创建测试VASP的ENCUT输入文件目录以及准备vasp的输入文件")
#    try:
    print("Note: --------------------")
    start, stop, step = sys.argv[1:]
    print("    python test_qe_ecutwfc.py 0.2 1.2 0.2")
    print("    python test_qe_ecutwfc.py 1.0 0.0 -0.2")
    print("    以上两种都是产生0.2到1.0， 间隔为0.2的kspacing")
    kresolutions = list(np.round(np.arange(eval(start), eval(stop), eval(step)), decimals=2))
#    except:
#        print("你设置的有问题，采用一下的k点密度")
#        kresolutions = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.1, 0.05]
    print("Note: --------------------")
    print("    该脚本不提供自动提任务的命令: 你可以用以下命令提供命令:")
    print("        for i in {}; do cd $i; qsub submit.sh;   cd ..; done".format(' '.join(list(map(str, kresolutions)))))
    print("        for i in {}; do cd $i; sbatch submit.sh; cd ..; done".format(' '.join(list(map(str, kresolutions)))))
    old_scf_in_path = os.path.abspath("scf.in")
    for kresolution in kresolutions:
        test_path = os.path.abspath(str(kresolution))
        new_scf_in_path = os.path.join(test_path, "scf.in")
        if not os.path.exists(test_path):
            os.mkdir(test_path)
        change_kpoints(old_scf_in_path, new_scf_in_path, kresolution)
        write_submit(jobsystem, test_path)


    print("    尝试获得每原子的焓值 Et(meV/atom)")
    kresolution_Et = []
    for kresolution in kresolutions:
        scf_out_path = os.path.join(str(kresolution), "scf.out")
        Et_per_atom = get_total_energy(scf_out_path)
        timetot = checktime(scf_out_path)
        kresolution_Et.append([kresolution, Et_per_atom, timetot])
    kresolution_Et = sorted(kresolution_Et, key=lambda x: x[0],  reverse=True)



    # 方式一：通过相邻两个kpoints的能量的差判断收敛性
    # Ediff = np.diff(kresolution_Et, axis=0)[:,-1]
    # Ediff = np.insert(Ediff, 0, [0.0])
    # Ediff = Ediff[:, np.newaxis]
    # kresolution_Et_Ediff = np.hstack((kresolution_Et, Ediff))

    # 方式二：所有kpoints的能量减去最密kpoints的能量判断收敛性
    with open("kresolution.csv", 'w') as f:
        print("{:<12},{:<14},{:<14},  {:<14},{:<14},  {}".format("kresolution", "Et(Ry/atom)", "diff(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)", "time"))
        f.write("{:<12},{:<14},{:<14},  {:<14},{:<14},  {}".format("kresolution", "Et(Ry/atom)", "diff(Ry/atom)", "Et(eV/atom)", "diff(meV/atom)", "time"))
        for idx, Et_per_atom in enumerate(kresolution_Et):
            kresolution = Et_per_atom[0]
            Et_per_at   = Et_per_atom[1]
            timetot     = Et_per_atom[2]
            Ediff       = Et_per_at - kresolution_Et[-1][1]
            f.write("{:<12.3f},{:<14.8f},{:<14.8f},  {:<14.8f},{:<14.8f},  {}\n".format(kresolution, Et_per_at, Ediff, Et_per_at*13.605693009, Ediff*13.605693009*1000, timetot))
            print("{:<12.3f},{:<14.8f},{:<14.8f},  {:<14.8f},{:<14.8f},  {}".format(kresolution, Et_per_at, Ediff, Et_per_at*13.605693009, Ediff*13.605693009*1000, timetot))
        else:
            print("All scf.outs are OK, kresolution.csv has been wroten in current position")
