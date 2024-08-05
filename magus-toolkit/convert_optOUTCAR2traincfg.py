#!/usr/bin/env python3

import os
import sys
import time
import fnmatch
import subprocess

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
        if "OUTCAR" in files :
            ourcar_path = system_name.joinpath("OUTCAR")
            # 抑制错误消息
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {ourcar_path}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
            #    success_d.append(system_name.__str__())
                success_d.append(os.path.abspath(system_name))
            else:
            #    fail_d.append(system_name.__str__())
                fail_d.append(os.path.abspath(system_name))
        elif "INCAR" in files or "POTCAR" in files or "POSCAR" in files:
            # none_d.append(system_name.__str__())
            none_d.append(os.path.abspath(system_name))


# def custom_sort(item):
#     parts = item.split("/")
#     return eval(parts[-1])

def custom_sort(item):
    parts = item.split("/")
    try:
        return (parts[:-1], eval(parts[-1]))
    except:
        return (parts[:-1], parts[-1])

print("\nrelax-succeeded")
success_d = sorted(success_d, key=custom_sort)
for line in success_d:
    print(line)


# 获取当前目录
current_dir = os.getcwd()

# 获取指定目录（如果有）并转换为绝对路径
if len(sys.argv) > 1:
    specified_dir = os.path.abspath(sys.argv[1])
else:
    specified_dir = ""  # 如果没有指定目录，则保持为空字符串

# 检查是否指定路径
if specified_dir:
    target_train_cfg = os.path.join(specified_dir, 'train.cfg')
    # 检查指定路径下的train.cfg是否存在
    if not os.path.exists(target_train_cfg):
        target_train_cfg = os.path.join(current_dir, 'train.cfg')
        print("train.cfg hasn't existed!")
        open(target_train_cfg, 'a').close()
    else:
        print(f"train.cfg has existed in {target_train_cfg}")
else:
    target_train_cfg = os.path.join(current_dir, 'train.cfg')
    # 检查当前路径下的train.cfg是否存在
    if not os.path.exists(target_train_cfg):
        target_train_cfg = os.path.join(current_dir, 'train.cfg')
        print("train.cfg hasn't existed!")
        open(target_train_cfg, 'a').close()
    else:
        print(f"train.cfg has existed in {target_train_cfg}")

# 遍历当前目录下的所有子目录
for root in success_d:
    # 切换到包含OUTCAR的目录
    os.chdir(root)
    # 执行命令
    command = ['mlp', 'convert-cfg', 'OUTCAR', 'train.cfg', '--input-format=vasp-outcar']
    subprocess.run(command)

    # 读取生成的train.cfg文件内容
    with open('train.cfg', 'r') as src_file:
        src_content = src_file.read()
    
    # 将内容追加到目标train.cfg文件中
    with open(target_train_cfg, 'a') as dest_file:
        dest_file.write(src_content)

# 切换回初始目录
os.chdir(current_dir)