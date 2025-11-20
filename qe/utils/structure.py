"""
结构处理工具模块

提供结构文件读取、对称性分析等功能。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Tuple, List, Dict
import numpy as np
from ase.io import read as ase_read
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from core.types import Matrix3x3, Vector3D, Composition


class StructureReader:
    """结构文件读取器"""

    @staticmethod
    def read_structure(file_path: Path) -> Structure:
        """
        读取结构文件

        支持格式：POSCAR, CONTCAR, CIF, XYZ 等

        Parameters
        ----------
        file_path : Path
            结构文件路径

        Returns
        -------
        Structure
            Pymatgen Structure对象
        """
        # 使用ASE读取，然后转换为Pymatgen
        atoms = ase_read(str(file_path))

        # 转换为Pymatgen Structure
        structure = Structure(
            lattice=atoms.get_cell(),
            species=atoms.get_chemical_symbols(),
            coords=atoms.get_scaled_positions(),
            coords_are_cartesian=False
        )

        return structure

    @staticmethod
    def get_symmetry_info(structure: Structure) -> Dict:
        """
        获取对称性信息

        Parameters
        ----------
        structure : Structure
            晶体结构

        Returns
        -------
        Dict
            对称性信息字典
        """
        sga = SpacegroupAnalyzer(structure)

        return {
            'space_group': sga.get_space_group_symbol(),
            'space_group_number': sga.get_space_group_number(),
            'point_group': sga.get_point_group_symbol(),
            'crystal_system': sga.get_crystal_system(),
        }

    @staticmethod
    def get_high_symmetry_kpath(structure: Structure) -> Dict:
        """
        获取高对称K点路径

        Parameters
        ----------
        structure : Structure
            晶体结构

        Returns
        -------
        Dict
            高对称点路径信息
        """
        from pymatgen.symmetry.bandstructure import HighSymmKpath

        kpath = HighSymmKpath(structure)

        return {
            'kpath': kpath.kpath,
            'kpoints': kpath.kpath['kpoints'],
            'path': kpath.kpath['path']
        }
