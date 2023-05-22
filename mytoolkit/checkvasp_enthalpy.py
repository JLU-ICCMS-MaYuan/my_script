#!/usr/bin/env python3

import os

for root, dirs, files in os.walk(os.getcwd()):
    if ("INCAR" in files) and ("POTCAR" in files) and ("POSCAR" in files):
        incar_path = os.path.join(root, "INCAR")
        ismear = int(os.popen(f"grep ISMEAR {incar_path}").read().strip('\n').split()[-1])
        print("{} ISMEAR={}".format(incar_path, ismear))
        if ismear == 1:
            os.system(f"sed -i 's/ISMEAR   = 1/ISMEAR   = 0/g' {incar_path}")