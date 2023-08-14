#!/usr/bin/env python3
import os
import time
import fnmatch

from pathlib import Path
from pprint import pprint


fail_d = []
success_d = []
none_d = []
for system_name in Path.cwd().iterdir():
    if system_name.is_dir():
        relax_out = system_name.joinpath("relax.out")
        if relax_out.exists():
            res = os.popen(f'grep -s "JOB DONE" {relax_out}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
                success_d.append(system_name.__str__())
            else:
                fail_d.append(system_name.__str__())
        else:
            none_d.append(system_name.__str__())
print("\nrelax-none")
for line in none_d:
    print(line)

print("\nrelax-succeeded")
for line in success_d:
    print(line)

print("\nrelax-failed")
for line in fail_d:
    print(line)


fail_d = []
success_d = []
none_d = []
for system_name in Path.cwd().iterdir():
    if system_name.is_dir():
        relax_out = system_name.joinpath("scf.out")
        if relax_out.exists():
            res = os.system(f'grep -s "JOB DONE" {relax_out}') # 如果没有找到指定内容不输出错误结果。
            if res != 0:
                fail_d.append(system_name.__str__())
            else:
                success_d.append(system_name.__str__())
        else:
            none_d.append(system_name.__str__())
print("\nscf-none")
for line in none_d:
    print(line)

print("\nscf-succeeded")
for line in success_d:
    print(line)

print("\nscf-failed")
for line in fail_d:
    print(line)



