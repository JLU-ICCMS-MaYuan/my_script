#!/usr/bin/env python3

import os
import sys

prefix = os.popen("grep 'prefix      =' epw_fermi_nest.in").read().strip().split('=')[1].replace("'","").replace(",", "")
print(prefix)

nest_fucntion = os.popen("grep 'Nesting function (q)= ' epw_fermi_nest.out | awk '{print $4}' ").readlines()
print("nest_number", len(nest_fucntion))

if os.path.exists(f"../{prefix}_band.dat"):
    high_symmetry_path = os.popen(f"head -n {len(nest_fucntion)} ../{prefix}_band.dat " + " | awk '{print $1}'").readlines()
elif os.path.exists(f"{prefix}_band.dat"):
    high_symmetry_path = os.popen(f"head -n {len(nest_fucntion)}    {prefix}_band.dat " + " | awk '{print $1}'").readlines()
else:
    print(f"No existance of {prefix}_band.dat. The script will exit.")
    sys.exit(1)
    
with open(f"{prefix}.nesting_fn.dat", "w") as f:
    f.write("high_symmetry_path  Nesting_function\n")
    for point, nest in zip(high_symmetry_path, nest_fucntion):
        f.write(point.strip() + '  ' + nest.strip() + '\n') 
print(f"It has been write in {prefix}.nesting_fn.dat")
print(f"The corresponding high symmetry name is in the {prefix}_band.labelinfo.dat")
