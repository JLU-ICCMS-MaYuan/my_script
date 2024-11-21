#!/usr/bin/env python3
import time
import os
import subprocess
import argparse
import math
import glob

header = """#!/bin/bash -l
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

submit = r"""
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

# 删除所有以 slurm_batch_group_ 开头的 .sh 文件
def delete_old_slurm_scripts():
    slurm_scripts = glob.glob('slurm_batch_group_*.sh')  # 查找所有匹配的文件
    for script in slurm_scripts:
        # print(f"Deleting old Slurm script: {script}")
        os.remove(script)

def get_groups(INPUT_FILE:str, MAX_RUNNING_JOBS:int) -> list[list[str]] :
    # 读取任务目录
    with open(INPUT_FILE, 'r') as file:
        directories = file.readlines()

    # 清理目录列表中的空行和换行符
    directories = [d.strip() for d in directories if d.strip()]

    # 如果没有任务目录，则返回空列表
    if not directories:
        print("Error: No directories found in the input file.")
        return []

    # 将任务分成MAX_RUNNING_JOBS组
    groups = [directories[i:i + MAX_RUNNING_JOBS] for i in range(0, len(directories), MAX_RUNNING_JOBS)]
    
    print(f"Each group contains {len(groups[0])} tasks")
    return groups

# 定义生成Slurm脚本的函数
def generate_slurm_script(group, group_index):
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

    args = parser.parse_args()
    MAX_RUNNING_JOBS = args.max_running_jobs  # 最大并行任务数
    INPUT_FILE = args.input_file  # 输入的任务目录文件

    # 将任务分成MAX_RUNNING_JOBS组
    groups = get_groups(INPUT_FILE, MAX_RUNNING_JOBS)

    # 删除旧的 Slurm 脚本
    delete_old_slurm_scripts()

    # 生成Slurm脚本并保存到文件中
    for index, group in enumerate(groups):
        slurm_script = generate_slurm_script(group, index)
        slurm_script_filename = f"slurm_batch_group_{index + 1}.sh"

        # 将脚本内容写入文件
        with open(slurm_script_filename, 'w') as slurm_file:
            slurm_file.write(slurm_script)

        # print(f"Generated Slurm script for job group {index + 1}: {slurm_script_filename}")

        # 提交作业
        submit_command = f"sbatch {slurm_script_filename}"
        subprocess.run(submit_command, shell=True)

    print(f"Total {len(groups)} job groups submitted.")
