#!/usr/bin/env python3
import time
import os
import subprocess
import argparse
import math
import glob

header_forChem = """#!/bin/bash -l
#SBATCH --job-name=vasp
#SBATCH --partition=phys
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=56
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --account=mslab

source /public/env/intel2021
ulimit -v unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SLURM_EXPORT_ENV=ALL

"""

header_forInspur = """#!/bin/bash -l
#SBATCH  --job-name=vasp
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
##SBATCH  --partition=intel6240r_384
#SBATCH  --partition=intel6240r_192
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         

source /work/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0

"""

header_forRiken = """#!/bin/sh                           
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

submit_forChem = r"""
for dir in "${dirs[@]}"; do
    echo "go into $dir"
    cd "$dir" || continue

    # 执行指定代码
    for i in {1..4}; do
        cp "INCAR_$i" "INCAR"
        killall -9 vasp_std
        sleep 3
        srun --mpi=pmi2 /public/software/vasp.6.3.0/bin/vasp_std > vasp.log_$i 2>&1
        cp CONTCAR POSCAR
    done
    killall -9 vasp_std
done
"""

submit_forInspur = r"""
for dir in "${dirs[@]}"; do
    echo "go into $dir"
    cd "$dir" || continue

    # 执行指定代码
    for i in {1..4}; do
        cp "INCAR_$i" "INCAR"
        killall -9 vasp_std
        sleep 3
        mpirun -np 48 /work/home/mayuan/software/vasp.6.1.0/bin/vasp_std > vasp.log_$i 2>&1
        cp CONTCAR POSCAR
    done
    killall -9 vasp_std
done
"""

submit_forRiken = r"""
for dir in "${dirs[@]}"; do
    echo "go into $dir"
    cd "$dir" || continue

    # 执行指定代码
    for i in {1..4}; do
        cp "INCAR_$i" "INCAR"
        killall -9 vasp_std
        sleep 3
        srun /lustre/home/h240012/soft/vasp.6.1.0/bin/vasp_std > vasp.log_$i 2>&1
        cp CONTCAR POSCAR
    done
    killall -9 vasp_std
done
"""
# 删除所有以 slurm_batch_group_ 开头的 .sh 文件
def delete_old_slurm_scripts():
    slurm_scripts = glob.glob('slurm_batch_group_*.sh')  # 查找所有匹配的文件
    for script in slurm_scripts:
        # print(f"Deleting old Slurm script: {script}")
        os.remove(script)

def get_groups(input_file:str, max_running_jobs:int) -> list[list[str]] :
    # 读取任务目录
    with open(input_file, 'r') as file:
        directories = file.readlines()

    # 清理目录列表中的空行和换行符
    directories = [d.strip() for d in directories if d.strip()]

    # 如果没有任务目录，则返回空列表
    if not directories:
        print("Error: No directories found in the input file.")
        return []

    # 将任务分成max_running_jobs组
    groups = [directories[i:i + max_running_jobs] for i in range(0, len(directories), max_running_jobs)]

    print(f"Each group contains {len(groups[0])} tasks")
    return groups

# 定义生成Slurm脚本的函数
def generate_slurm_script(group, group_index, task_type):
    if task_type == "inspur":
        header = header_forInspur
        submit = submit_forInspur
    elif task_type == "chem":
        header = header_forChem
        submit = submit_forChem
    elif task_type == "riken":
        header = header_forRiken
        submit = submit_forRiken

    slurm_script = header
    slurm_script += f"# Job group {group_index + 1}\n"
    slurm_script += "echo \"Starting the group job...\"\n"
    dirs_paths = "dirs=(\n"  + '\n'.join(group) + "\n)"
    slurm_script = slurm_script + dirs_paths + submit
    # 返回生成的脚本内容
    return slurm_script

if __name__ == "__main__":
    print("You can use by: python batchsubmit.py -i relax-failed -m 10")

    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Submit jobs with controlled parallelism.")
    parser.add_argument('-m', "--max_running_jobs", type=int, default=10, help="Maximum number of jobs to run in parallel.")
    parser.add_argument('-i', "--input_file", type=str, help="File containing the directories to process.", required=True)
    parser.add_argument('-t', "--task_type", type=str, help="Task type and machine type. Supported: inspur, chem, riken", required=True)
    parser.add_argument('-s', "--submit", action='store_true', help="If set, submit the jobs; otherwise, just generate the scripts.")

    args = parser.parse_args()
    max_running_jobs = args.max_running_jobs  # 最大并行任务数
    input_file = args.input_file  # 输入的任务目录文件
    task_type  = args.task_type
    submit = args.submit

    # 将任务分成max_running_jobs组
    groups = get_groups(input_file, max_running_jobs)

    # 删除旧的 Slurm 脚本
    delete_old_slurm_scripts()

    # 生成Slurm脚本并保存到文件中
    for index, group in enumerate(groups):
        slurm_script = generate_slurm_script(group, index, task_type)
        slurm_script_filename = f"slurm_batch_group_{index + 1}.sh"

        # 将脚本内容写入文件
        with open(slurm_script_filename, 'w') as slurm_file:
            slurm_file.write(slurm_script)

        # print(f"Generated Slurm script for job group {index + 1}: {slurm_script_filename}")

        # 如果设置了 --submit 参数，提交作业
        if submit:
            submit_command = f"sbatch {slurm_script_filename}"
            subprocess.run(submit_command, shell=True)

    print(f"Total {len(groups)} job groups submitted.")
