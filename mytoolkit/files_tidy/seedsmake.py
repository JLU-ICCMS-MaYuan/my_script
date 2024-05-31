#!/usr/bin/env python3
import os
import re

from pathlib import Path

def makeseeds(seeds_dir, poscars_dir):

    if Path(seeds_dir).exists():
        with open(seeds_dir, 'r') as f:
            lines = f.readlines()
            ids = [i for i, line in enumerate(lines) if re.search(r'POSCAR-0\d{3}', line)]
    else:
        ids = []
        print(f"Attention! {seeds_dir} doesn't exist!")

    if not Path(poscars_dir).exists():
        raise FileExistsError(f"{poscars_dir} doesn't exist!")

    for posid, file in enumerate(Path(poscars_dir).glob(r"POSCAR*")):
        with open(file, 'r') as pos:
            poscar = pos.readlines()
        poscar[0] = "POSCAR-{:0>{}}\n".format(posid+1+len(ids), 4)
        with open(seeds_dir, 'a') as seed:
            seed.writelines(poscar)

if __name__ == "__main__":
    seeds_dir = input("指定一个种子文件\n")
    poscars_dir = input("指定你准备的POSCAR文件目录\n")
    makeseeds(seeds_dir, poscars_dir)
