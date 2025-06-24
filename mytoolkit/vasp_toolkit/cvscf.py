#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

fail_d = []
success_d = []
none_d = []

def parse_args():
    parser = argparse.ArgumentParser(description="Check SCF convergence with directory depth limit.")
    parser.add_argument('-d', '--depth', type=int, default=None,
                        help='Specify max depth for directory search (default: unlimited)')
    return parser.parse_args()

def walk_with_depth(base_dir, max_depth):
    base_dir = Path(base_dir).resolve()
    for root, dirs, files in os.walk(base_dir):
        rel_path = Path(root).relative_to(base_dir)
        if max_depth is not None and len(rel_path.parts) > max_depth:
            # Prune deeper directories
            dirs[:] = []
            continue
        yield root, dirs, files

def custom_sort(item):
    parts = item.split("/")
    try:
        return (parts[:-1], eval(parts[-1]))
    except:
        return parts

def main():
    args = parse_args()

    for root, dirs, files in walk_with_depth(".", args.depth):
        system_name = Path(root)
        if "OUTCAR" in files:
            outcar_path = system_name / "OUTCAR"
            converge1 = os.popen(f'grep -s -a "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" {outcar_path}').read()
            converge2 = os.popen(f'grep -s -a "General timing and accounting informations for this job" {outcar_path}').read()
            kbstar    = os.popen(f'grep -s -a "in kB  \*"  {outcar_path}').read()
            if converge1 and converge2 and not kbstar:
                nelm = os.popen(f'grep -s -a "NELM   =" {outcar_path} | cut -d \';\' -f 1 | cut -d \'=\' -f 2').read()
                iteration = os.popen(f'grep -s -a "Iteration" {outcar_path} | wc -l').read()
                try:
                    if int(iteration.strip()) < int(nelm.strip()):
                        success_d.append(str(system_name.resolve()))
                    else:
                        fail_d.append(str(system_name.resolve()))
                except:
                    fail_d.append(str(system_name.resolve()))
            else:
                fail_d.append(str(system_name.resolve()))
        else:
            none_d.append(str(system_name.resolve()))

    # 输出和保存文件
    for label, data, filename in [("scf-none", none_d, "scf-none"),
                                  ("scf-succeeded", success_d, "scf-succeeded"),
                                  ("scf-failed", fail_d, "scf-failed")]:
        print(f"\n{label}")
        try:
            data = sorted(data, key=custom_sort)
        except:
            pass
        for line in data:
            print(line)
        with open(filename, 'w') as f:
            f.write('\n'.join(data))

if __name__ == "__main__":
    main()
