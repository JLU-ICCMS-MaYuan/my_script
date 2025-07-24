#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phonon band structure using Phonopy API

支持软件: VASP / QE / YAML
依赖: phonopy >=3.0, numpy, matplotlib (可选)

使用说明：
1. 设置 soft（vasp/qe/yaml）
2. 设置结构文件 / 力常数文件 / 超胞维度
3. 设置高对称路径点和 labels
"""

from typing import Sequence
import numpy as np
import phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.interface.qe import PH_Q2R, read_pwscf
from phonopy.phonon.band_structure import get_band_qpoints
import matplotlib.pyplot as plt

# ------------------ 载入 phonon ------------------ #

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
    return phonon

# ------------------ 运行 band 结构 ------------------ #

def run_band_structure(phonon, path_points, labels):

    phonon.run_band_structure(path_points, with_eigenvectors=False, labels=labels)
    band_plot = phonon.plot_band_structure()
    band_plot.savefig("phonon_band.png", dpi=300)

# ------------------ 主函数入口 ------------------ #

if __name__ == "__main__":
    soft = "qe"   # "vasp", "qe", or "yaml"
    structure_file = "scf.in"
    fc_file = "Ce1Sc2H24.fc"
    supercell_matrix = [6, 6, 6]

    # 高对称路径设置（含每段插值点数）
    path_points = [
        [0.0, 0.0, 0.0],         # Gamma
        [0.5, 0.0, 0.0],         # M
        [1/3, 1/3, 0.0],         # K
        [0.0, 0.0, 0.0],         # Gamma
        [0.0, 0.0, 0.5],         # A
        [0.5, 0.0, 0.5],         # L
        [1/3, 1/3, 0.5],         # H
        [0.0, 0.0, 0.5],         # A
        [0.5, 0.0, 0.5],         # L
        [0.5, 0.0, 0.0],         # M
        [1/3, 1/3, 0.0],         # K
        [1/3, 1/3, 0.5],         # H
    ]
    labels = ["Γ", "M", "K", "Γ", "A", "L", "H", "A", "L", "M", "K", "H"]
    points = get_band_qpoints([path_points], 101)
    phonon = load_phonon(soft, structure_file, fc_file, supercell_matrix)
    distances, frequencies = run_band_structure(phonon, points, labels)

