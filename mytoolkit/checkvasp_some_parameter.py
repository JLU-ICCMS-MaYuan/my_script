#!/usr/bin/env python3

import os

for root, dirs, files in os.walk(os.getcwd()):
    if ("INCAR" in files) and ("POTCAR" in files) and ("POSCAR" in files):
        incar_path = os.path.join(root, "INCAR")
        ismear = int(os.popen(f"grep NCORE {incar_path}").read().strip('\n').split()[-1])
        print("{} NCORE={}".format(incar_path, ismear))
        if ismear == 1:
            os.system(f"sed -i 's/NCORE    = 1/NCORE    = 4/g' {incar_path}")