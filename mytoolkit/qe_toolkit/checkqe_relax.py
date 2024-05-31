#!/usr/bin/env python3
import os
import time
import fnmatch

from pathlib import Path
from pprint import pprint

fail_d = []
success_d = []
without_relaxout = []
print("Note: --------------------")
print("    检查qe批量结构优化是否完成, 讲结构优化成果的任务写入当前目录的relax-succeeded.dat文件中")
print("    默认每个结构优化的压强为 300GPa, 是否改变, 如果改变请输入你的目标压强(输入空格或者回车意味着采用默认压强300GPa")

pressure = input() or "300" # 给 input 函数设置默认值.如果用户没有输入任何内容，就会使用默认值

for system_name in Path.cwd().iterdir():
    if system_name.is_dir():
        work_pressure = system_name.joinpath(pressure)
        relax_out = work_pressure.joinpath("relax.out")
        if relax_out.exists():
            res = os.system(f'grep -s "JOB DONE" {relax_out}') # 如果没有找到指定内容不输出错误结果。
            if res != 0:
                fail_d.append(work_pressure.__str__())
            else:
                success_d.append(work_pressure.__str__())
        else:
            without_relaxout.append(work_pressure.__str__())

with open("relax-without.dat", "w") as f:
    for line in without_relaxout:
        f.write(line+"\n")

with open("relax-succeeded.dat", "w") as f:
    for line in success_d:
        f.write(line+"\n")

with open("relax-failed.dat", "w") as f:
    for line in fail_d:
        f.write(line+"\n")





