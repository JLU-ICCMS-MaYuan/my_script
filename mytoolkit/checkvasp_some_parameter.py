#!/usr/bin/env python3

import os

for root, dirs, files in os.walk(os.getcwd()):
    if ("INCAR" in files) and ("POTCAR" in files) and ("POSCAR" in files):
        incar_path = os.path.join(root, "INCAR")
        encut = int(os.popen(f"grep ENCUT {incar_path}").read().strip('\n').split()[-1])
        if encut == 800:
            print(root)