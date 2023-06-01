#!/usr/bin/env python3

import os

for root, dirs, files in os.walk(os.getcwd()):
    if ("INCAR" in files) and ("POTCAR" in files) and ("POSCAR" in files):
        if "OUTCAR" in files:
            line = os.popen("grep 'reached required accuracy - stopping structural energy minimisation' {}".format(os.path.join(root, "OUTCAR"))).read().split()
            if line:
                print(f"{root} is ok")
            else:
                print(f"{root} is something wrong")
        else:
            print(f"{root} doesn't has OUTCAR")
