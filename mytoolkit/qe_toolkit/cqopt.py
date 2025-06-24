#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Check job status (relax.out or scf.out) with directory depth limit.")
    parser.add_argument("-d", "--depth", type=int, default=None, help="Max directory depth to search (default: unlimited)")
    return parser.parse_args()

def walk_with_depth(base_dir, max_depth):
    base_dir = Path(base_dir).resolve()
    for root, dirs, files in os.walk(base_dir):
        rel_path = Path(root).relative_to(base_dir)
        if max_depth is not None and len(rel_path.parts) > max_depth:
            dirs[:] = []  # Stop deeper search
            continue
        yield Path(root), files

def check_job_status(files, root_path, filename):
    file_path = root_path / filename
    if filename in files:
        res = os.popen(f'grep -s "JOB DONE" {file_path}').read()
        if res.strip():
            return "succeeded"
        else:
            return "failed"
    return "none"

def categorize_by_file(filename, depth):
    result = {
        "succeeded": [],
        "failed": [],
        "none": []
    }

    for root_path, files in walk_with_depth(".", depth):
        status = check_job_status(files, root_path, filename)
        result[status].append(str(root_path.resolve()))

    return result

def write_and_print(label, data):
    print(f"\n{label}")
    for line in data:
        print(line)
    with open(label, 'w') as f:
        f.write('\n'.join(data))

def main():
    args = parse_args()

    # Check relax.out
    relax_result = categorize_by_file("relax.out", args.depth)
    write_and_print("relax-none", relax_result["none"])
    write_and_print("relax-succeeded", relax_result["succeeded"])
    write_and_print("relax-failed", relax_result["failed"])

    # Check scf.out
    scf_result = categorize_by_file("scf.out", args.depth)
    write_and_print("scf-none", scf_result["none"])
    write_and_print("scf-succeeded", scf_result["succeeded"])
    write_and_print("scf-failed", scf_result["failed"])

if __name__ == "__main__":
    main()
