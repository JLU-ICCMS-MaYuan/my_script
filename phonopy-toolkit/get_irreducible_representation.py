#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Irreducible Representation Analysis for Arbitrary q-points using Phonopy API

使用说明
--------
1. 将该脚本放置于包含 `phonopy_params.yaml` 或相关结构和力常数文件的目录中。
2. 根据所用软件（VASP/QE/YAML）修改 `soft` 变量，并设置目标 q 点。
3. 执行脚本获取该 q 点下每个模式对应的不可约表示和频率信息。

依赖：phonopy >=3.0, spglib, sympy
"""

from pathlib import Path
from typing import Sequence

import numpy as np
import phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.interface.qe import PH_Q2R, read_pwscf
from phonopy.phonon.irreps import IrReps


def load_phonon(soft: str, structure_file: str, fc_file: str, supercell_matrix: Sequence[int]):
    if soft == "vasp":
        cell = read_vasp(structure_file)
        phonon = phonopy.load(
            supercell_matrix=supercell_matrix,
            unitcell=cell,
            force_sets_filename=fc_file,
            calculator="vasp",
        )
    elif soft == "qe":
        cell, _ = read_pwscf(structure_file)
        q2r = PH_Q2R(fc_file)
        q2r.run(cell)
        q2r.write_force_constants()
        phonon = phonopy.load(
            supercell_matrix=supercell_matrix,
            unitcell=cell,
            force_constants_filename="force_constants.hdf5",
            calculator="qe",
        )
    elif soft == "yaml":
        phonon = phonopy.load(structure_file)
    else:
        raise ValueError("soft must be one of 'vasp', 'qe', or 'yaml'")

    print("[INFO] Phonopy object loaded.")
    return phonon


def analyze_irreps(phonon, qpt: Sequence[float]):
    print(f"[INFO] Analyzing irreducible representations at q-point: {qpt}")
    irreps = IrReps(phonon.dynamical_matrix, qpt, symprec=1e-5, degeneracy_tolerance=1e-5)
    irreps.run()
    frequencies = phonon.get_frequencies(qpt)
    labels = irreps._get_ir_labels()  # 推荐替换 _get_ir_labels() 为属性方式
    deg_sets = irreps._degenerate_sets
    print(labels)
    print(deg_sets)
    
    print("\n[INFO] Irreducible Representations (with degeneracies):")
    for band_idx, deg_group in enumerate(deg_sets):
        # 所有简并模式的编号、频率、ir label
        modes_idx = []; freqs = []
        for i in deg_group:
            modes_idx.append(f"Mode {i+1}")
            freqs.append(frequencies[i])
        ir_label = labels[band_idx]   # 简并组共用一个 ir_label

        print(
            f"{' & '.join(modes_idx):20s} | "
            f"{ir_label:<6s} | "
            f"Freq: {freqs[0]:7.3f} THz | {freqs[0]*33.35641:7.2f} cm⁻¹ | "
            f"Degeneracy: {len(deg_group)}"
        )



if __name__ == "__main__":
    # ---------- 用户输入区 ---------- #
    soft = "qe"  # 选择 'vasp', 'qe' 或 'yaml'
    structure_file = "scf.in"          # 对于 VASP 则为 POSCAR-init；对于 YAML 则为 phonopy_params.yaml
    fc_file = "Ce1Sc2H24.fc"          # 对于 VASP 为 FORCE_SETS；对于 YAML 可不填
    supercell_matrix = [6, 6, 6]       # 超胞维度

    qpt = [  0,   0,   0]  # Γ 点
    # qpt = [1/2, 0.0, 0.0]  # M 点
    # qpt = [0.0, 0.0, 1/2]  # A 点
    # qpt = [1/3, 1/3,   0]  # K 点
    # qpt = [1/2,   0, 1/2]  # L 点

    # ---------- 执行分析 ---------- #
    phonon = load_phonon(soft, structure_file, fc_file, supercell_matrix)
    analyze_irreps(phonon, qpt)

