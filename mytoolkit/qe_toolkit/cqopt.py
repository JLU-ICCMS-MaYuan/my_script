#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

fail_d = []
success_d = []
none_d = []

def parse_args():
    parser = argparse.ArgumentParser(description="Check relaxation status with directory depth limit.")
    parser.add_argument('-d', '--depth', type=int, default=None,
                        help='Specify max depth for directory search (default: unlimited)')
    return parser.parse_args()

def walk_with_depth(base_dir, max_depth):
    base_dir = Path(base_dir).resolve()
    for root, dirs, files in os.walk(base_dir):
        rel_path = Path(root).relative_to(base_dir)
        if max_depth is not None and len(rel_path.parts) > max_depth:
            dirs[:] = []  # Prune deeper traversal
            continue
        yield root, dirs, files

def custom_sort(item):
    parts = item.split("/")
    try:
        return (parts[:-1], str(eval(parts[-1])))
    except:
        return (parts[:-1], parts[-1])

def main():
    args = parse_args()

    for root, dirs, files in walk_with_depth(".", args.depth):
        system_name = Path(root)
        files_set = set(files)
        if "OUTCAR" in files_set:
            outcar_path = system_name / "OUTCAR"
            res = os.popen(f'grep -s "reached required accuracy - stopping structural energy minimisation" {outcar_path}').read()
            if res.strip():
                success_d.append(str(system_name.resolve()))
            else:
                fail_d.append(str(system_name.resolve()))
        elif {"INCAR", "POTCAR", "POSCAR"} & files_set:
            none_d.append(str(system_name.resolve()))

    for label, data, filename in [("relax-none", none_d, "relax-none"),
                                  ("relax-succeeded", success_d, "relax-succeeded"),
                                  ("relax-failed", fail_d, "relax-failed")]:
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
