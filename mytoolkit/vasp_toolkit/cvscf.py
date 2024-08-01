#!/usr/bin/env python3

import os
import time
import fnmatch

from pathlib import Path
from pprint import pprint


fail_d = []
success_d = []
none_d = []
# for system_name in Path.cwd().iterdir():
    # if system_name.is_dir():
        # ourcar_path = system_name.joinpath("OUTCAR")
for root, dirs, files in os.walk("."):
        system_name = Path(root)
        if "OUTCAR" in files:
            ourcar_path = system_name.joinpath("OUTCAR")
            # 抑制错误消息
            converge1 = os.popen(f'grep -s -a "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" {ourcar_path}').read()
            converge2 = os.popen(f'grep -s -a "General timing and accounting informations for this job" {ourcar_path}').read()
            kbstar   = os.popen(f'grep -s -a "in kB  \*"  {ourcar_path}').read()
            if converge1 and converge2 and not kbstar:
                nelm = os.popen(f'grep -s -a "NELM   =" {ourcar_path} ' + " | cut -d ';' -f 1 | cut -d '=' -f 2 ").read()
                iteration = os.popen(f'grep -s -a "Iteration" {ourcar_path} | wc -l').read() # 如果没有找到指定内容不输出错误结果。
                try:
                    if int(iteration) < int(nelm):
                        # success_d.append(system_name.__str__())
                        success_d.append(os.path.abspath(system_name))
                    else:
                        # fail_d.append(system_name.__str__())
                        fail_d.append(os.path.abspath(system_name))
                except:
                    # fail_d.append(system_name.__str__())
                    fail_d.append(os.path.abspath(system_name))
            else:
                # fail_d.append(system_name.__str__())
                fail_d.append(os.path.abspath(system_name))
        else:
            # none_d.append(system_name.__str__())
            none_d.append(os.path.abspath(system_name))


# def custom_sort(item):
#     parts = item.split("/")
#     return eval(parts[-1])

def custom_sort(item):
    parts = item.split("/")
    return (parts[:-1], eval(parts[-1]))

print("\nscf-none")
for line in none_d:
    print(line)

print("\nscf-succeeded")
try: 
    success_d = sorted(success_d, key=custom_sort)
except:
    success_d = success_d
for line in success_d:
    print(line)

print("\nscf-failed")
try: 
    fail_d = sorted(fail_d, key=custom_sort)
except:
    fail_d = fail_d
for line in fail_d:
    print(line)
