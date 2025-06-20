#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phonon modulation workflow (VASP / QE).

使用方法
--------
1. 把本脚本放到与计算产生的输入文件同一目录。
   - VASP: 需要 `POSCAR-init`、`FORCE_SETS`
   - QE  : 需要 `scf.in`、`*.fc` (此处默认 `Ce1Sc2H24.fc`)
   - YAML: 需要 `phonopy_params.yaml`
2. 根据需要修改 `soft`（'vasp' / 'qe' / 'yaml'）以及入口处的调制参数。
3. 运行： `python modulation.py`
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import os
import phonopy
import numpy as np

from phonopy.interface.vasp import read_vasp, write_vasp
from phonopy.interface.qe import PH_Q2R, read_pwscf
from phonopy.interface.vasp import read_vasp
from phonopy.interface.calculator import get_default_cell_filename
# --------------- 构建 Phonopy ---------------- #


def load_phonon_from_vasp(
    poscar_in: str = "POSCAR-init",
    force_sets: str = "FORCE_SETS",
    supercell_matrix: list[int] | None = None,
):
    cell = read_vasp(poscar_in)

    # 2. 处理力常数
    fc_file = Path(force_sets)

    # 3. 构建 Phonopy 对象
    phonon = phonopy.load(
        supercell_matrix=supercell_matrix,
        unitcell=cell,
        force_sets_filename=fc_file,
        calculator="vasp",
    )

    # 4. 保存 yaml 备份
    phonon.save(settings={"force_constants": True})
    print("[INFO] Phonopy object built from VASP files.")
    return phonon


def load_phonon_from_qe(
    scf_in: str = "scf.in",
    fc_raw: str = "Ce1Sc2H24.fc",
    supercell_matrix: Sequence[int] | None = None,
):
    """
    从 QE 输出读取并返回 Phonopy 对象。

    Parameters
    ----------
    scf_in : str
        QE `pw.x` 的输入文件（包含晶胞信息）
    fc_raw : str
        QE `matdyn.x` 或 `q2r.x` 生成的力常数文件
    supercell_matrix : list[int]
        超胞维度
    """
    cell, _ = read_pwscf(scf_in)

    # 将 .fc 转换为 Phonopy 可读的 force_constants.hdf5
    q2r = PH_Q2R(fc_raw)
    q2r.run(cell)
    q2r.write_force_constants()
    phonon = phonopy.load(
        supercell_matrix=supercell_matrix,
        calculator="qe",
        unitcell=cell,
        force_constants_filename="force_constants.hdf5",
    )
    # 保存 yaml 备份，今后可直接读取
    phonon.save(settings={"force_constants": True})
    print("[INFO] Phonopy object built from QE files.")
    return phonon


def load_phonon_from_yaml(yaml_file: str = "phonopy_params.yaml"):
    """
    直接从 phonopy YAML 恢复 Phonopy 对象。

    Parameters
    ----------
    yaml_file : str
        `phonon.save()` 生成的 YAML 文件
    """
    print(f"[INFO] Loading Phonopy object from {yaml_file}.")
    return phonopy.load(yaml_file)

# --------------- 调制计算 ---------------- #


def run_modulation(
    soft,
    phonon,
    qpt: Sequence[float] | np.ndarray,
    amplitude_list: Sequence[float] | np.ndarray,
    nmode: int,
    dim: Sequence[int],
    natoms_in_primitive: int | None = None,
    
) -> None:
    """
    计算并写出调制后的 POSCAR & YAML。

    Parameters
    ----------
    phonon : Phonopy
        已具备力常数的 Phonopy 实例
    qpt : length‑3
        q 点（分数坐标）
    amplitude_list : iterable
        调制振幅
    nmode : int
        带指数（0 表示第一条带）
    dim : length‑3
        超胞维度
    natoms_in_primitive : int | None
        原胞原子数；缺省取 phonon.unitcell 大小
    """
    Na = natoms_in_primitive or len(phonon.unitcell)

    # Phonopy 的 band index 从 0 开始
    phonon_modes = [
        [list(qpt), nmode, amp * np.sqrt(Na), 0.0] for amp in amplitude_list
    ]

    phonon.run_modulations(dim, phonon_modes)
    phonon.write_yaml_modulations()
    # phonon.write_modulations()
    
    base_fname = get_default_cell_filename(interface_mode=None)
    deltas = []
    for amp, u in zip(amplitude_list, phonon._modulation._u):
        cell = phonon._modulation._get_cell_with_modulation(u)
        write_vasp(filename=f"M{base_fname}_{(amp):2.2f}", cell=cell)
        deltas.append(u)
    
        sum_of_deltas = np.sum(deltas, axis=0)
        cell = phonon._modulation._get_cell_with_modulation(sum_of_deltas)
        write_vasp(filename=f"M{base_fname}", cell=cell)
        
        no_modulations = np.zeros(sum_of_deltas.shape, dtype=complex)
        cell = phonon._modulation._get_cell_with_modulation(no_modulations)
        write_vasp(filename=f"M{base_fname}-orig", cell=cell)

        
    if soft == "qe":
        os.system('sed -i "s/   1.0/    0.529177/g" MPOSCAR*')
    print(f"[INFO] Modulation finished → {Path.cwd()}")
    

# --------------- 主函数 ---------------- #


def main(
    qpt: Sequence[float],
    amplitude_list: Sequence[float],
    nmode: int,
    dim: Sequence[int],
    input_structure: str,
    input_fc_name: str,
    input_supercell_matrix: Sequence[float],
    soft,
) -> None:
    """
    指定调制参数，然后调用 run_modulation。
    若以后接入 argparse，只需把四个形参换成解析结果即可。
    """
    if soft == "vasp":
        phonon = load_phonon_from_vasp(
            poscar_in=input_structure,
            force_sets=input_fc_name,
            supercell_matrix=input_supercell_matrix,
        )
    elif soft == "qe":
        phonon = load_phonon_from_qe(
            scf_in=input_structure,
            fc_raw=input_fc_name,
            supercell_matrix=input_supercell_matrix,
        )
    elif soft == "yaml":
        phonon = load_phonon_from_yaml()
    else:
        raise ValueError("soft must be 'vasp', 'qe' or 'yaml'.")

    run_modulation(
        soft=soft,
        phonon=phonon,
        qpt=list(qpt),
        amplitude_list=np.asarray(amplitude_list),
        nmode=nmode,
        dim=list(dim),
    )


# --------------- 脚本入口 ---------------- #

if __name__ == "__main__":
    # 在此处手动设置调制参数；如需命令行解析，替换成 argparse 即可
    qpt_   = [0., 0., 0.]          # Γ 点
    amps_  = np.linspace(start=-10, stop=10, num=21)      # 三个振幅
    nmode_ = 0                     # 0 代表的是第一条带，最底下的那条带，1代表的是第二条带，这是因为python是从0开始计数
    dim_   = [1, 1, 1]                # 2×2×2 超胞
    #input_structure_ = "POSCAR-init"
    #input_fc_name_ = "FORCE_SETS"
    #input_supercell_matrix_ = [2., 2., 2.]
    #soft = 'vasp'
    input_structure_ = "scf.in"
    input_fc_name_ = "Ce1Sc2H24.fc"
    input_supercell_matrix_ = [3., 3., 3.]
    soft='qe'
    #soft='yaml'
    main(
        qpt=qpt_, 
        amplitude_list=amps_, 
        nmode=nmode_, 
        dim=dim_, 
        input_structure=input_structure_, 
        input_fc_name=input_fc_name_, 
        input_supercell_matrix=input_supercell_matrix_,
        soft=soft,
        )
