"""
电子性质绘图模块

提供电子能带、态密度、投影能带的高质量绘图。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Optional, Dict
import numpy as np
import matplotlib.pyplot as plt

from analysis.plotting.base import BasePlotter


class ElectronBandPlotter(BasePlotter):
    """电子能带绘图器"""

    def plot(self, kpoints: np.ndarray, energies: np.ndarray,
             fermi_energy: float = 0.0,
             projected_weights: Optional[Dict] = None,
             output_path: Optional[Path] = None,
             **kwargs) -> Path:
        """
        绘制电子能带（支持投影）

        Parameters
        ----------
        kpoints : np.ndarray
            K点数组
        energies : np.ndarray
            能量数组 (nkpts, nbands)
        fermi_energy : float
            费米能级 [eV]
        projected_weights : Dict, optional
            投影权重 {元素+轨道: weights数组}
        output_path : Path

        Returns
        -------
        Path
        """
        fig, ax = plt.subplots(figsize=(8, 6))

        # 调整能量到费米能级
        energies_shifted = energies - fermi_energy

        if projected_weights is None:
            # 普通能带图
            for i in range(energies.shape[1]):
                ax.plot(kpoints, energies_shifted[:, i], 'b-', linewidth=1.0)
        else:
            # 投影能带图（散点图，点大小表示权重）
            for orb, weights in projected_weights.items():
                for i in range(energies.shape[1]):
                    # 归一化权重到合适的点大小
                    sizes = weights[:, i] * 100
                    ax.scatter(kpoints, energies_shifted[:, i],
                             s=sizes, alpha=0.6, label=orb if i == 0 else "")

        # 费米能级
        ax.axhline(y=0, color='k', linestyle='--', linewidth=1.5, label='$E_F$')

        # 标签
        ax.set_xlabel('Wave vector', fontsize=14)
        ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
        ax.set_title('Electronic Band Structure', fontsize=16, pad=15)

        if projected_weights:
            ax.legend(loc='best', frameon=False, fontsize=10)

        ax.grid(True, alpha=0.3, linestyle=':')

        # 保存
        if output_path:
            self.save_figure(fig, output_path)
            plt.close()
            return output_path

        plt.show()
        return None


class ElectronDOSPlotter(BasePlotter):
    """电子态密度绘图器"""

    def plot(self, energies: np.ndarray, dos: np.ndarray,
             fermi_energy: float = 0.0,
             projected_dos: Optional[Dict] = None,
             output_path: Optional[Path] = None,
             **kwargs) -> Path:
        """
        绘制电子态密度（支持投影）

        Parameters
        ----------
        energies : np.ndarray
            能量数组
        dos : np.ndarray
            总态密度
        fermi_energy : float
            费米能级
        projected_dos : Dict, optional
            投影态密度 {元素+轨道: PDOS数组}
        output_path : Path

        Returns
        -------
        Path
        """
        fig, ax = plt.subplots(figsize=(6, 6))

        # 调整能量
        energies_shifted = energies - fermi_energy

        # 总态密度
        ax.plot(dos, energies_shifted, 'k-', linewidth=2.0, label='Total DOS')

        # 投影态密度（堆叠或填充）
        if projected_dos:
            colors = plt.cm.tab10(np.linspace(0, 1, len(projected_dos)))
            for (label, pdos), color in zip(projected_dos.items(), colors):
                ax.fill_betweenx(energies_shifted, 0, pdos,
                                alpha=0.5, label=label, color=color)

        # 费米能级
        ax.axhline(y=0, color='k', linestyle='--', linewidth=1.5)

        # 标签
        ax.set_xlabel('DOS (states/eV)', fontsize=14)
        ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
        ax.set_title('Electronic Density of States', fontsize=16, pad=15)

        if projected_dos:
            ax.legend(loc='best', frameon=False, fontsize=10)

        ax.grid(True, alpha=0.3, linestyle=':')

        # 保存
        if output_path:
            self.save_figure(fig, output_path)
            plt.close()
            return output_path

        plt.show()
        return None
