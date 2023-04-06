#!/usr/bin/env python3

import re
import os
import sys
import shutil
from pprint import pprint

from pathlib import Path

def split_seeds(filename, dst_dir):
    '''
    Args:
        filenanme: the name of seed file
        dst_dir: the position where all seeds will be stored.
    '''
    if not Path(dst_dir).exists():
        os.mkdir(dst_dir)
    if not Path(filename).exists():
        raise FileExistsError(f"{filename} doesn't exist!")
    with open(filename) as f:
        lines = f.readlines()
        ids = [i for i, line in enumerate(lines) if re.search(r'POSCAR-0\d{3}', line)]
        ids.append(len(lines))
        names_paths = []
        for i in range(len(ids)-1):
            i1 = ids[i]
            i2 = ids[i+1]
            struct_lst = lines[i1:i2]
            seed_file_path = Path(dst_dir).joinpath(struct_lst[0].strip('\n')+'.vasp')
            with open(seed_file_path, 'w') as poscar:
                poscar.writelines(struct_lst)
            names_paths.append([struct_lst[0].strip('\n'), seed_file_path])

    return names_paths

def prepare_opt(press, pre_dirpath, names_paths):
    if not Path(pre_dirpath).exists():
        raise FileExistsError(f"{pre_dirpath} doesn't exist!")
    cwd = Path.cwd() 
    for name, seed_file_path in names_paths:
        print(f"prepare {name}")
        press_name_dstpath = cwd.joinpath(str(press), name)
        if not press_name_dstpath.exists():
            os.makedirs(press_name_dstpath)
        slurm_path = Path(pre_dirpath).joinpath('slurm.sh')
        poscar_path = Path(press_name_dstpath).joinpath("POSCAR")
        shutil.copy(seed_file_path, press_name_dstpath)
        shutil.copy(slurm_path, press_name_dstpath)
        shutil.copy(seed_file_path, press_name_dstpath)
        shutil.copy(seed_file_path, poscar_path)
        for i in range(1,5):
            incar_srcpath = Path(pre_dirpath).joinpath(f'INCAR_{str(i)}')
            incar_dstpath = Path(press_name_dstpath).joinpath(f'INCAR_{str(i)}')
            shutil.copy(incar_srcpath, press_name_dstpath)
            incar_dstpath = str(incar_dstpath)
            os.system(f"sed -i 's/PSTRESS = 2000/PSTRESS = {float(press)*10}/g' {str(incar_dstpath)}")

        eles=os.popen(f"sed -n '6p' {poscar_path}").read().split()
        potcar_dstfile= Path(press_name_dstpath).joinpath("POTCAR")
        if potcar_dstfile.exists():
            os.remove(potcar_dstfile)
        for e in eles:
            potcar_srcdir = Path.cwd().joinpath("potcar_lib", e)
            os.system(f"cat {potcar_srcdir} >> {potcar_dstfile}")

    return "Done"


if __name__ == '__main__':
    print("你必须制定四个输入参数")
    fn      = input("种子文件的路径\n").strip()
    dd      = input("拆分后种子文件的存储位置(你准备将种子文件中的结构文件拿出来以后都放到哪里?)\n").strip()
    press   = input("结构优化压强(GPa, 指定该目录路径后会在该路径下生成一个press命名的目录, 所有的结构优化都在该目录下进行)\n").strip()
    pre_dir = input("VASP输入文件的准备目录的路径(一个名叫prepare的目录, 里面包含了4个INCAR_*和一个slurm.sh)\n").strip()

    # fn      = "seeds.vasp"
    # dd      = "struct_lib"
    # press   = "200.0"
    # pre_dir = "prepare"

    names_paths = split_seeds(fn, dd)
    result = prepare_opt(press, pre_dir, names_paths)
    print(result)

