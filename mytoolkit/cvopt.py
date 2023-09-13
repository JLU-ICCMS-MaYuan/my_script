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
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {ourcar_path}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
                success_d.append(system_name.__str__())
            else:
                fail_d.append(system_name.__str__())
        else:
            none_d.append(system_name.__str__())


# def custom_sort(item):
#     parts = item.split("/")
#     return eval(parts[-1])

def custom_sort(item):
    parts = item.split("/")
    try:
        return (parts[:-1], eval(parts[-1]))
    except:
        return (parts[:-1], parts[-1])

print("\nrelax-none")
for line in none_d:
    print(line)

print("\nrelax-succeeded")
success_d = sorted(success_d, key=custom_sort)
for line in success_d:
    print(line)

print("\nrelax-failed")
fail_d = sorted(fail_d, key=custom_sort)
for line in fail_d:
    print(line)





