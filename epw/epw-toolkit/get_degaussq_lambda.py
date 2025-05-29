#!/usr/bin/env python3

import os
import re

def get_prefix_and_grid():
    prefix = None
    nkf = [None, None, None]
    nqf = [None, None, None]

    for fname in ["epw_iso_sc.in", "epw_aniso_sc.in"]:
        if os.path.exists(fname):
            with open(fname, 'r') as f:
                for line in f:
                    # 匹配 prefix
                    match_prefix = re.match(r"\s*prefix\s*=\s*['\"]?([^'\"\s]+)", line)
                    if match_prefix:
                        prefix = match_prefix.group(1)

                    # 匹配 nkf1, nkf2, nkf3
                    for i in range(1, 4):
                        match_nkf = re.match(rf"\s*nkf{i}\s*=\s*(\d+)", line)
                        if match_nkf:
                            nkf[i - 1] = match_nkf.group(1)

                    # 匹配 nqf1, nqf2, nqf3
                    for i in range(1, 4):
                        match_nqf = re.match(rf"\s*nqf{i}\s*=\s*(\d+)", line)
                        if match_nqf:
                            nqf[i - 1] = match_nqf.group(1)
            break

    if not prefix:
        raise ValueError("未在输入文件中找到 prefix")
    if None in nkf or None in nqf:
        raise ValueError("未能在输入文件中找到 nkf1~3 或 nqf1~3 参数")

    return prefix, nkf, nqf

def extract_elph_data(input_file, output_file, nkf, nqf):
    coupling_values = []
    smearing_values = []

    with open(input_file, 'r') as f:
        lines = f.readlines()

    read_coupling = False
    read_smearing = False

    for line in lines:
        line = line.strip()

        if line.startswith("Integrated el-ph coupling"):
            read_coupling = True
            continue
        if read_coupling and line.startswith("#"):
            coupling_values = line.split()[1:]
            read_coupling = False
            continue

        if line.startswith("Phonon smearing"):
            read_smearing = True
            continue
        if read_smearing and line.startswith("#"):
            smearing_values = line.split()[1:]
            read_smearing = False
            continue

    if len(coupling_values) != len(smearing_values):
        raise ValueError("Smearing 和 Coupling 数量不一致！")

    with open(output_file, 'w') as f:
        title = f"λ (nkf={' '.join(nkf)}, nqf={' '.join(nqf)})"
        f.write(f"Phonon_smearing(meV)\t{title}\n")
        for s, c in zip(smearing_values, coupling_values):
            f.write(f"{s}\t{c}\n")


if __name__ == "__main__":
    prefix, nkf, nqf = get_prefix_and_grid()
    input_filename = f"{prefix}.a2f"
    output_filename = f"{prefix}_degaussq_lambda.dat"
    extract_elph_data(input_filename, output_filename, nkf, nqf)
    print(f"✅ 已生成文件: {output_filename}")

