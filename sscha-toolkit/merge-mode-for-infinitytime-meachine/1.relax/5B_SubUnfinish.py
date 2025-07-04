from numpy import *
import sys
import subprocess
import numpy as np
import time

import os
def get_que_num():
    """获取当前目录下正在运行的任务数量"""
    current_directory = os.getcwd()  # 获取当前工作目录
    
    # 使用 squeue 命令并过滤包含当前目录的任务
    command = f"squeue -h --format '%Z' | grep {current_directory}"
    try:
        squeue_output = subprocess.check_output(command, shell=True).decode('utf-8').strip()
        # print(squeue_output)  # 打印当前队列输出，便于调试
        # 返回匹配的行数（即正在运行的任务数量）
        return len(squeue_output.splitlines()), squeue_output.splitlines()
    except subprocess.CalledProcessError:
        #print('Error fetching queue data')  # 错误处理
        return 0, []  # 如果 squeue 没有输出或其他错误，返回 0

num = int(sys.argv[1])
directory = "run_calculation"
output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly

print("You get number of output",len(output_files))
print("Check these output right or wrong")
energies = np.zeros(len(output_files))
id_nums = []
unfinish = []
for file in output_files:
    # Get the number of the configuration.
    id_number = int(file.split("_")[-1].split(".")[0])
    id_nums.append(id_number)

    # Load the file
    ff = open(file, "r")
    lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
    ff.close()

    Flag_Ener=False
    for l in lines:
        if len(l) > 0 :
           if l.split()[0] == "!":
              Flag_Ener=True
              Flag_stress=False
              for l in lines:
                  if len(l) > 0 :
                     if l.split()[0] == "total" and l.split()[1] == "stress":
                        Flag_stress=True

              if Flag_stress==False :
                 print("Stress WRONG",id_number)
                 unfinish.append(id_number)


    if Flag_Ener==False :
       unfinish.append(id_number)
       print("ENERGY WRONG",id_number)


sort_num = sorted(id_nums)
sort_num.append(10897654)
j=0
for i in range(0,num):
    if i==sort_num[j]:
       j=j+1
    else:
       print("without",i)
       unfinish.append(i)

# print(unfinish)
n_sub = len(unfinish)

#all_scf_files = [os.path.join("run_calculation", f) for f in os.listdir("run_calculation") if f.startswith("espresso_run_")]
header="""#!/bin/sh 
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=intel6430
#SBATCH  --nodes=1                          
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                         

source /data/home/mayuan/intel-2024/oneapi/intel2024_setvars.sh --force

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ulimit -s unlimited
"""

MAX_RUNNING_JOBS = 4
#for id_unfinish in unfinish:
while True:
    que_num, queue_path = get_que_num()
    print(que_num)
    # 从任务队列中取出下一个任务目录
    if que_num < MAX_RUNNING_JOBS:
        id_unfinish = unfinish.pop(0)
        filename= "un_sub_{}.sh".format(i)
        run_line="mpirun -np 64 /data/home/mayuan/soft/qe-7.3.1-intel-oneAPI2024/bin/pw.x  -nk 4 -in  ESP_{}.pwi > ESP_{}.pwo 2>&1".format(id_unfinish,id_unfinish)
        with open(filename, "w") as f:
            f.write(header)
            print(run_line,file=f)
        os.system(rf"mv {filename} run_calculation && cd run_calculation && sbatch {filename} && cd ../")
    else:
        time.sleep(20)
