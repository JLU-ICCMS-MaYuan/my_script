#!/usr/bin/env python3

import os
import re
import sys
import glob
import argparse
import subprocess

from pathlib import Path
from pprint import pprint

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', default='', help='please input where you want to save total_train.cfg.')
parser.add_argument('-e', '--extract-last-configuration', action="store_true", help='whether or not you only want to save the last configurations')
args = parser.parse_args()

fail_d = []
success_d = []
none_d = []
# for system_name in Path.cwd().iterdir():
    # if system_name.is_dir():
        # ourcar_path = system_name.joinpath("OUTCAR")
for root, dirs, files in os.walk("."):
    system_name = Path(root)
    outcar_files = glob.glob(str(system_name.joinpath("OUTCAR*")))
    if outcar_files:
        for outcar_path in outcar_files:
            # 抑制错误消息
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {outcar_path}').read() # 如果没有找到指定内容不输出错误结果。
            if res:
                print(outcar_path)
                success_d.append(os.path.abspath(outcar_path))
            else:
                fail_d.append(os.path.abspath(outcar_path))
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

def extract_specified_configuration(extract_last_configuration):
    with open('train.cfg', 'r') as src_file:
        src_content = src_file.read()

    if extract_last_configuration:
        #print(src_content); input()
        pattern = r'BEGIN_CFG.*?END_CFG'
        matches = re.findall(pattern, src_content, re.DOTALL)
        # 如果有匹配项，提取最后一个匹配项的内容
        try:
            last_match = matches[-1]
            return last_match
        except:
            return src_content
    else:
        return src_content

print("\nrelax-succeeded")
success_d = sorted(success_d, key=custom_sort)
for line in success_d:
    print(line)
with open('relax-succeeded', 'w') as f:
    f.writelines('\n'.join(success_d))

# 获取当前目录
current_dir = os.getcwd()

# 获取指定目录（如果有）并转换为绝对路径
if args.path:
    total_train_cfg = os.path.join(args.path, 'total_train.cfg')
else:
    total_train_cfg = os.path.join(current_dir, 'total_train.cfg')

if os.path.exists(total_train_cfg):
    os.system(f"rm -fr {total_train_cfg}")

# 遍历当前目录下的所有子目录
for succ in success_d:
    # 执行命令
    command = ['mlp', 'convert-cfg', succ, 'train.cfg', '--input-format=vasp-outcar']
    subprocess.run(command)
    src_content = extract_specified_configuration(args.extract_last_configuration)
    # 将内容追加到目标train.cfg文件中
    with open(total_train_cfg, 'a') as dest_file:
        dest_file.write(src_content)

# 切换回初始目录
os.chdir(current_dir)
