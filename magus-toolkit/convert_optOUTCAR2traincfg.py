#!/usr/bin/env python3

import os
import re
import sys
import glob
import argparse
import subprocess

from pathlib import Path

def custom_sort(item):
    parts = item.split("/")
    try:
        return (parts[:-1], eval(parts[-1]))
    except:
        return (parts[:-1], parts[-1])

def get_all_succeeded_OUTCAR(pressures):
    fail_d = []
    none_d = []
    success_d = []
    for root, dirs, files in os.walk("."):
        system_name = Path(root)
        outcar_files = glob.glob(str(system_name.joinpath("OUTCAR*")))
        if outcar_files:
            for outcar_path in outcar_files:
                res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {outcar_path}').read()
                if res:
                    pressure = str(int(float(os.popen(f'grep -s "PSTRESS="  {outcar_path} |  ' + "awk '{print $2}'").read().strip()) / 10))
                    if pressure in pressures:
                        print(outcar_path)
                        success_d.append(os.path.abspath(outcar_path))
                else:
                    fail_d.append(os.path.abspath(outcar_path))
        elif "INCAR" in files or "POTCAR" in files or "POSCAR" in files:
            # none_d.append(system_name.__str__())
            none_d.append(os.path.abspath(system_name))

    print("\n relaxed-succeeded")
    success_d = sorted(success_d, key=custom_sort)
    for line in success_d:
        print(line)
    with open('relaxed-succeeded', 'w') as f:
        f.writelines('\n'.join(success_d))
    return success_d

# 废弃的函数，不是因为写粗了，而是因为找到了mlp获得多个构型最后一帧的命令行参数--last
def extract_specified_configuration(extract_last_configuration):
    with open('train.cfg', 'r') as src_file:
        src_content = src_file.read()

    if extract_last_configuration:
        #print(src_content); input()
        pattern = r'BEGIN_CFG.*?END_CFG'
        matches = re.findall(pattern, src_content, re.DOTALL)
        # 如果有匹配项，提取最后一个匹配项的内容
        try:
            last_match = matches[-1] + '\n'
            return last_match
        except:
            return src_content
    else:
        return src_content

def write_train_cfg(success_d, total_train_cfg):
    if os.path.exists(total_train_cfg):
        os.system(f"rm -fr {total_train_cfg}")

    # 遍历当前目录下的所有子目录
    for succ in success_d:
        # 执行命令
        if args.extract_last_configuration:
            command = ['mlp', 'convert-cfg', succ.strip('\n'), 'train.cfg', '--input-format=vasp-outcar', '--last']
        else:
            command = ['mlp', 'convert-cfg', succ.strip('\n'), 'train.cfg', '--input-format=vasp-outcar']
        print('  '.join(command))
        subprocess.run(command)
        # src_content = extract_specified_configuration(args.extract_last_configuration)
        # 将内容追加到目标train.cfg文件中
        with open('train.cfg', 'r') as src_file:
            src_content = src_file.read()
        with open(total_train_cfg, 'a') as dest_file:
            dest_file.write(src_content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--filename',  required=True, default='total_train.cfg', help='please input where you want to save total_train.cfg.')
    parser.add_argument('-p', '--pressures',  required=True, default=[], nargs='+', type=str, help='please input the pressure, which is named as filename')
    parser.add_argument('-l', '--extract-last-configuration',action="store_true", help='whether or not you only want to save the last configurations')

    args = parser.parse_args()

    # 获取指定目录（如果有）并转换为绝对路径
    current_dir = os.getcwd()
    
    total_train_cfg = os.path.join(args.filename)
    
    # 获取当前目录
    if os.path.exists("relaxed-succeeded"):
        os.system("rm relaxed-succeeded")
        success_d = get_all_succeeded_OUTCAR(args.pressures)
    else:
        success_d = get_all_succeeded_OUTCAR(args.pressures)

    write_train_cfg(success_d, total_train_cfg)
    # 切换回初始目录
    os.chdir(current_dir)

