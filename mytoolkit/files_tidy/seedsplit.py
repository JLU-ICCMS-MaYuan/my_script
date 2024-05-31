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
    if not Path.cwd().joinpath(dst_dir).exists():
        os.mkdir(dst_dir)
    if not Path.cwd().joinpath(filename).exists():
        print("stlib doesn't exist!")
        sys.exit(1)

    with open(filename) as f:
        lines = f.readlines()
        ids = [i for i, line in enumerate(lines) if re.search(r":", line)]
        # ids = [i for i, line in enumerate(lines) if re.search("[0-9]+-[a-zA-Z]+-*", line)]
        ids.append(len(lines))
        names_paths = []
        for i in range(len(ids)-1):
            i1 = ids[i]
            i2 = ids[i+1]
            struct_lst = lines[i1:i2]
            print("title is ", struct_lst[0])
            seed_file_path = Path(dst_dir).joinpath(str(i+1)+"-Y-Al-H"+'.vasp')
            with open(seed_file_path, 'w') as poscar:
                poscar.writelines(struct_lst)
            names_paths.append([struct_lst[0].strip('\n'), seed_file_path])

    return names_paths


if __name__ == '__main__':

    print("Note: --------------------")
    print("    注意在使用该脚本是, 请务必确保你的种子文件的名称为seeds.vasp, 然后所有的种子会被拆分为XXX.vasp。所以请注意种子文件中的结构格式为POSCAR格式。")
    names_paths = split_seeds("seeds.vasp", "stlib")