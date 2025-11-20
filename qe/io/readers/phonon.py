"""
声子输出文件读取器

读取QE声子计算的输出文件。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import re


class PhononOutputReader:
    """声子输出文件读取器"""

    @staticmethod
    def read_frequencies(freq_file: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取声子频率文件

        Parameters
        ----------
        freq_file : Path
            频率文件路径（如 matdyn.freq）

        Returns
        -------
        qpoints : np.ndarray
            Q点坐标数组
        frequencies : np.ndarray
            频率数组 (nqpts, nbands)
        """
        with open(freq_file, 'r') as f:
            lines = f.readlines()

        # 解析第一行获取维度
        match = re.search(r'nbnd=\s*(\d+),\s*nks=\s*(\d+)', lines[0])
        if match:
            nbands = int(match.group(1))
            nqpts = int(match.group(2))
        else:
            raise ValueError(f"无法解析频率文件头部: {lines[0]}")

        qpoints = []
        frequencies = []

        i = 1
        while i < len(lines):
            # 读取Q点坐标
            parts = lines[i].strip().split()
            if len(parts) == 3:
                qpoint = [float(x) for x in parts]
                qpoints.append(qpoint)

                # 读取该Q点的所有频率
                i += 1
                freqs = []
                while i < len(lines) and len(lines[i].strip().split()) > 3:
                    freq_vals = [float(x) for x in lines[i].strip().split()]
                    freqs.extend(freq_vals)
                    i += 1

                frequencies.append(freqs[:nbands])
            else:
                i += 1

        return np.array(qpoints), np.array(frequencies)

    @staticmethod
    def read_phonon_dos(dos_file: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取声子态密度文件

        Parameters
        ----------
        dos_file : Path
            态密度文件路径

        Returns
        -------
        frequencies : np.ndarray
            频率数组
        dos : np.ndarray
            态密度数组
        """
        data = np.loadtxt(dos_file, skiprows=1)
        frequencies = data[:, 0]  # 第一列是频率
        dos = data[:, 1]          # 第二列是DOS

        return frequencies, dos

    @staticmethod
    def read_gamma_lines(gamma_file: Path) -> Dict[float, np.ndarray]:
        """
        读取声子线宽文件 (gam.lines)

        Parameters
        ----------
        gamma_file : Path
            Gamma线宽文件路径

        Returns
        -------
        Dict[float, np.ndarray]
            {展宽值: 线宽数组}
        """
        gamma_data = {}

        with open(gamma_file, 'r') as f:
            lines = f.readlines()

        current_broaden = None
        current_gamma = []

        for line in lines:
            # 检查是否是展宽值行
            if 'Broadening' in line:
                # 保存上一个展宽的数据
                if current_broaden is not None and current_gamma:
                    gamma_data[current_broaden] = np.array(current_gamma)
                    current_gamma = []

                # 读取新的展宽值
                match = re.search(r'([\d.]+)', line)
                if match:
                    current_broaden = float(match.group(1))

            # 读取数据行
            elif line.strip() and not line.strip().startswith('Gamma'):
                try:
                    values = [float(x) for x in line.strip().split()]
                    current_gamma.extend(values)
                except ValueError:
                    pass

        # 保存最后一个展宽的数据
        if current_broaden is not None and current_gamma:
            gamma_data[current_broaden] = np.array(current_gamma)

        return gamma_data
