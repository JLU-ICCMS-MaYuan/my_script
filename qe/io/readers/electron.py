"""
电子性质输出文件读取器

读取QE电子能带、态密度等输出文件。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Tuple, Optional
import numpy as np
import re


class ElectronOutputReader:
    """电子输出文件读取器"""

    @staticmethod
    def read_fermi_energy(scf_out: Path) -> float:
        """
        从SCF输出文件读取费米能级

        Parameters
        ----------
        scf_out : Path
            SCF输出文件路径

        Returns
        -------
        float
            费米能级 [eV]
        """
        with open(scf_out, 'r') as f:
            content = f.read()

        # 搜索费米能级
        match = re.search(r'the Fermi energy is\s+([\d.]+)\s+ev', content, re.IGNORECASE)
        if match:
            return float(match.group(1))

        # 如果没找到，搜索最高占据态
        match = re.search(r'highest occupied level.*?:\s+([\d.]+)', content)
        if match:
            return float(match.group(1))

        raise ValueError(f"无法从{scf_out}读取费米能级")

    @staticmethod
    def read_band_structure(bands_out: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取能带结构数据

        Parameters
        ----------
        bands_out : Path
            能带输出文件路径

        Returns
        -------
        kpoints : np.ndarray
            K点数组
        energies : np.ndarray
            能量数组 (nkpts, nbands)
        """
        with open(bands_out, 'r') as f:
            lines = f.readlines()

        # 解析文件
        # TODO: 实现完整的能带文件解析
        return np.array([]), np.array([])

    @staticmethod
    def read_dos(dos_file: Path) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        读取态密度文件

        Parameters
        ----------
        dos_file : Path
            DOS文件路径

        Returns
        -------
        energies : np.ndarray
            能量数组
        dos : np.ndarray
            态密度数组
        fermi_energy : float
            费米能级
        """
        with open(dos_file, 'r') as f:
            # 读取第一行获取费米能级
            first_line = f.readline()
            match = re.search(r'EFermi\s*=\s*([\d.]+)', first_line)
            fermi_energy = float(match.group(1)) if match else 0.0

            # 读取数据
            f.seek(0)
            data = np.loadtxt(f, skiprows=1)

        energies = data[:, 0]  # 能量
        dos = data[:, 1]       # DOS

        return energies, dos, fermi_energy
