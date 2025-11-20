"""
输入文件写入器基类

提供QE输入文件生成的通用方法，消除代码重复。

作者：Claude
创建时间：2025-11-19
"""

from typing import List, Dict, TextIO
from pathlib import Path
from abc import ABC, abstractmethod

from core.types import Vector3D, Matrix3x3


class BaseWriter(ABC):
    """QE输入文件写入器基类"""

    @abstractmethod
    def write(self, output_path: Path) -> Path:
        """
        生成输入文件

        Parameters
        ----------
        output_path : Path
            输出文件路径

        Returns
        -------
        Path
            生成的文件路径
        """
        pass

    def write_atomic_species(self, f: TextIO, composition: Dict[str, int],
                            pseudopotentials: Dict[str, str]):
        """
        写入ATOMIC_SPECIES部分

        Parameters
        ----------
        f : TextIO
            文件对象
        composition : Dict[str, int]
            元素组成 {元素: 数量}
        pseudopotentials : Dict[str, str]
            赝势文件 {元素: 赝势文件名}
        """
        f.write("ATOMIC_SPECIES\n")
        for element in composition.keys():
            mass = self._get_atomic_mass(element)
            pp_file = pseudopotentials.get(element, f"{element}.UPF")
            f.write(f"  {element:5} {mass:10.5f} {pp_file}\n")
        f.write("\n")

    def write_cell_parameters(self, f: TextIO, lattice: Matrix3x3, unit: str = "angstrom"):
        """
        写入CELL_PARAMETERS部分

        Parameters
        ----------
        f : TextIO
            文件对象
        lattice : Matrix3x3
            晶格向量矩阵
        unit : str
            单位 ('angstrom', 'bohr', 'alat')
        """
        f.write(f"CELL_PARAMETERS {{{unit}}}\n")
        for row in lattice:
            f.write(f"  {row[0]:16.10f} {row[1]:16.10f} {row[2]:16.10f}\n")
        f.write("\n")

    def write_atomic_positions(self, f: TextIO, species: List[str],
                              positions: List[Vector3D], coord_type: str = "crystal"):
        """
        写入ATOMIC_POSITIONS部分

        Parameters
        ----------
        f : TextIO
            文件对象
        species : List[str]
            原子种类列表
        positions : List[Vector3D]
            原子坐标列表
        coord_type : str
            坐标类型 ('crystal', 'angstrom', 'bohr')
        """
        f.write(f"ATOMIC_POSITIONS {{{coord_type}}}\n")
        for sp, pos in zip(species, positions):
            f.write(f"  {sp:5} {pos[0]:16.10f} {pos[1]:16.10f} {pos[2]:16.10f}\n")
        f.write("\n")

    def write_k_points(self, f: TextIO, kgrid: tuple, shift: tuple = (0, 0, 0)):
        """
        写入K_POINTS部分（自动网格）

        Parameters
        ----------
        f : TextIO
            文件对象
        kgrid : tuple
            K点网格 (nk1, nk2, nk3)
        shift : tuple
            K点偏移
        """
        f.write("K_POINTS {automatic}\n")
        f.write(f"  {kgrid[0]} {kgrid[1]} {kgrid[2]}  {shift[0]} {shift[1]} {shift[2]}\n")
        f.write("\n")

    @staticmethod
    def _get_atomic_mass(element: str) -> float:
        """获取原子质量（简化版）"""
        from pymatgen.core import Element
        return float(Element(element).atomic_mass)
