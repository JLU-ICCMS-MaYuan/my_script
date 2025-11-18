#!/usr/bin/env python3
import re
import os
import sys

cwd = os.getcwd()
squeuelines = os.popen(r"squeue -u liuhanyu -h --format '%A    %t    %R    %M    %Z '"  + f" | grep {cwd}").readlines()

with open('log', 'r') as f:
    loglines = f.read()
matches1 = re.findall(r'job: ([a-f0-9]+) submit; job_id is (\d+)', loglines)
matches2 = re.findall(r'job:([a-f0-9]+) re-submit after terminated; new job_id is (\d+)', loglines)

hash_jobid_map = {match[1]: match[0] for match in matches1 + matches2}
# 正则表达式  
  
# 使用正则表达式搜索日志  
running = []
waiting = []
for sqline in squeuelines:
    #print(sqline)
    sqline = sqline.strip().split()
    jobid = sqline[0]
    time  = sqline[3]
    fppath = sqline[-1]
    hashid = hash_jobid_map[jobid]
    hash_subrun_path = os.path.join(fppath, hashid+'.sub.run')
    with open(hash_subrun_path, 'r') as f:
        hash_subrun_content = f.read()
    try:
        taskid = re.search(r'task\.\d+\.\d+', hash_subrun_content).group()
        task_path = os.path.join(fppath, taskid)
        if 'R' == sqline[1]:
            running.append([jobid, hashid, time, task_path])
        elif 'PD' == sqline[1]:
            waiting.append([jobid, hashid, time, task_path])
    except AttributeError:
        taskid = re.search(r'task\.\d+', hash_subrun_content).group()
        task_path = os.path.join(fppath, taskid)
        if 'R' == sqline[1]:
            running.append([jobid, hashid, time, task_path])
        elif 'PD' == sqline[1]:
            waiting.append([jobid, hashid, time, task_path])
    except Exception as e:
        pass



print('Waiting')
for waitid, hashid, time, waitpath in waiting:
    print(waitid, hashid, time, waitpath)
print('Running')
for runid, hashid, time, runpath in running:
    print(runid, hashid, time, runpath)