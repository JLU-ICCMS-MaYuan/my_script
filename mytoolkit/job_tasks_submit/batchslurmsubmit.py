#!/usr/bin/env python3
import time
import os
import sys
import subprocess
import argparse

print("You can use by: python batchsubmit.py -i relax-failed -m 10 -s 'jiaoben.sh'")

# 设置命令行参数解析
parser = argparse.ArgumentParser(description="Submit jobs with controlled parallelism.")
parser.add_argument('-m', "--max_running_jobs", type=int, help="Maximum number of jobs to run in parallel.")
parser.add_argument('-s', "--submit_command", type=str, help="Command to submit each job (e.g.,jiaoben.sh).")
parser.add_argument('-i', "--input_file", type=str, help="File containing the directories to process.")

args = parser.parse_args()

# 获取命令行参数
MAX_RUNNING_JOBS = args.max_running_jobs  # 最大并行任务数
SUBMIT_COMMAND = args.submit_command      # 提交任务的脚本
INPUT_FILE = args.input_file              # 任务目录列表文件

# 获取当前工作目录
cwd = os.getcwd()

# 读取任务目录列表
with open(INPUT_FILE, 'r') as f:
    start_que = [i.strip() for i in f.readlines()]  # 读取任务目录列表

# 当前正在运行的任务队列
run_que = []

def sub_job(pth):
    """提交任务的函数"""
    os.system(rf'cp {SUBMIT_COMMAND} {pth}')
    os.system(rf'cd {pth} && sbatch {SUBMIT_COMMAND} && cd {cwd}')

def get_que_num():
    """获取当前目录下正在运行的任务数量"""
    current_directory = os.getcwd()  # 获取当前工作目录
    
    # 使用 squeue 命令并过滤包含当前目录的任务
    command = f"squeue -h --format '%Z' | grep {current_directory}"
    try:
        squeue_output = subprocess.check_output(command, shell=True).decode('utf-8').strip()
        # print(squeue_output)  # 打印当前队列输出，便于调试
        # 返回匹配的行数（即正在运行的任务数量）
        print(len(squeue_output.splitlines()), squeue_output.splitlines())
        return len(squeue_output.splitlines()), squeue_output.splitlines()
    except subprocess.CalledProcessError:
        #print('Error fetching queue data')  # 错误处理
        return 0, []  # 如果 squeue 没有输出或其他错误，返回 0

if len(start_que) == 0:
    sys.exit(1)

# 检查是否还可以提交新的任务
while True:
    # 获取当前队列中的作业数目
    que_num, queue_path = get_que_num()

    if que_num < MAX_RUNNING_JOBS:

        # 从任务队列中取出下一个任务目录
        next_dir = start_que[0]
        
        # 检查路径是否已经在当前的队列中
        if next_dir in queue_path:
            # print(f"Skipping {next_dir} as it is already in the queue.")
            start_que.pop(0)  # 从队列中移除该路径，继续检查下一个路径
            continue
        
        # 将任务目录添加到运行队列
        run_que.append(start_que.pop(0))
        pth = run_que[-1]

        # 提交任务
        # print(f"Submitting job for directory: {pth}")
        sub_job(pth)
        
    if start_que == 0:
        break

    # 每20秒检查一次
    time.sleep(20)
