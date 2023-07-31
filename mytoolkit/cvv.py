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
        ourcar_path = system_name.joinpath("OUTCAR")
        if ourcar_path.exists():
            # 抑制错误消息
            res = os.system(f'grep -s "reached required accuracy - stopping structural energy minimisation" {ourcar_path}') # 如果没有找到指定内容不输出错误结果。
            if res != 0:
                fail_d.append(system_name.__str__())
            else:
                success_d.append(system_name.__str__())
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





