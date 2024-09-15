#!/usr/bin/env python3

import os
import sys
import re
import subprocess

def get_free_energy(outcar_path):
    try:
        result = subprocess.run("grep -s " + " free  energy   TOTEN " + f" {outcar_path} | tail -n 1 | awk '{{print $5}}'", 
                                shell=True, capture_output=True, text=True)
        dE = result.stdout.strip()
        return float(dE) if dE else 1e14
    except Exception:
        return 1e14

def get_num_atoms(outcar_path):
    try:
        result = subprocess.run(f'grep -n -s "position of ions in cartesian coordinates" {outcar_path}', 
                                shell=True, capture_output=True, text=True)
        begin_id = int(result.stdout.split(":")[0])

        N = 0
        with open(outcar_path, 'r') as f:
            for i, line in enumerate(f, start=1):
                if i > begin_id:
                    coords = re.findall(r"[-+]?\d+\.\d+", line.strip())
                    if len(coords) == 3:
                        N += 1
                    else:
                        break
        return N
    except Exception:
        return 1

def main():
    try:
        outcar_path = sys.argv[1]
    except IndexError:
        outcar_path = "OUTCAR"
    
    dE = get_free_energy(outcar_path)
    N = get_num_atoms(outcar_path)
    
    print("{:<12.8f} {:<3}".format(dE/N, N))

if __name__ == "__main__":
    main()
