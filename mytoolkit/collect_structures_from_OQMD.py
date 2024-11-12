#!/usr/bin/env python
import sys
from tqdm import tqdm
from pathlib import Path

from ase.io import read, write


def isequal_species(species1:list[str], species2:list[str]):
    """判断species1中全部元素是否都在species2中"""
    if len(species1) == len(species2):
        for spe in species1:
            if spe not in species2:
                return False
        else:
            return True

species = sys.argv[1:]

# 获取当前目录下的所有文件
source_dir = Path('oqmd_configs')
dstinate_dir = Path('_'.join(species))
dstinate_dir.mkdir(parents=True, exist_ok=True)
#if not dstinate_dir.exists():
#    os.mkdir(dstinate_dir)


# 筛选出以 .vasp 结尾的文件名
vasp_files = [file for file  in source_dir.iterdir() if file.suffix == '.vasp']

i = 0
for vf in tqdm(vasp_files, desc="Processing files"):
    at = read(vf)
    if isequal_species(species1=species, species2=set(at.symbols)):
        i += 1
        newname = dstinate_dir.joinpath(str(i)+'.'+f'{at.symbols}' + '.vasp')
        write(filename=newname, images=at, format='vasp')
