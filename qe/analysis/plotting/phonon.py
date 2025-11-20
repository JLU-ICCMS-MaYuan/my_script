"""
声子谱绘图模块

提供声子能带、声子态密度的高质量绘图功能。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Optional, List, Dict
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analysis.plotting.base import BasePlotter


class PhononBandPlotter(BasePlotter):
    """声子能带绘图器"""

    def plot(self, qpoints: np.ndarray, frequencies: np.ndarray,
             high_sym_points: Optional[Dict] = None,
             output_path: Optional[Path] = None,
             show_gamma: bool = True,
             **kwargs) -> Path:
        """
        绘制声子能带图

        Parameters
        ----------
        qpoints : np.ndarray
            Q点坐标数组 (npoints,)
        frequencies : np.ndarray
            频率数组 (npoints, nbands)
        high_sym_points : Dict, optional
            高对称点信息 {position: label}
        output_path : Path, optional
            输出文件路径
        show_gamma : bool
            是否显示Gamma线宽

        Returns
        -------
        Path
            输出文件路径
        """
        fig, ax = plt.subplots(figsize=(8, 6))

        # 绘制每条能带
        nbands = frequencies.shape[1]
        for i in range(nbands):
            ax.plot(qpoints, frequencies[:, i], 'b-', linewidth=1.0)

        # 设置高对称点
        if high_sym_points:
            for pos, label in high_sym_points.items():
                ax.axvline(x=pos, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
                ax.text(pos, ax.get_ylim()[1]*0.95, label,
                       ha='center', va='top', fontsize=14)

        # 零线
        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

        # 标签
        ax.set_xlabel('Wave vector', fontsize=14)
        ax.set_ylabel('Frequency (cm$^{-1}$)', fontsize=14)
        ax.set_title('Phonon Band Structure', fontsize=16, pad=15)

        # 网格
        ax.grid(True, alpha=0.3, linestyle=':')

        # 保存
        if output_path:
            self.save_figure(fig, output_path)
            plt.close()
            return output_path

        plt.show()
        return None


class PhononDOSPlotter(BasePlotter):
    """声子态密度绘图器"""

    def plot(self, frequencies: np.ndarray, dos: np.ndarray,
             projected_dos: Optional[Dict[str, np.ndarray]] = None,
             output_path: Optional[Path] = None,
             **kwargs) -> Path:
        """
        绘制声子态密度图

        Parameters
        ----------
        frequencies : np.ndarray
            频率数组
        dos : np.ndarray
            总态密度
        projected_dos : Dict, optional
            投影态密度 {元素: DOS数组}
        output_path : Path

        Returns
        -------
        Path
            输出文件路径
        """
        fig, ax = plt.subplots(figsize=(6, 6))

        # 总态密度
        ax.plot(dos, frequencies, 'k-', linewidth=2.0, label='Total')

        # 投影态密度
        if projected_dos:
            colors = plt.cm.tab10(np.linspace(0, 1, len(projected_dos)))
            for (element, pdos), color in zip(projected_dos.items(), colors):
                ax.plot(pdos, frequencies, linewidth=1.5, label=element, color=color)

        # 零线
        ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)

        # 标签
        ax.set_xlabel('Phonon DOS (states/cm$^{-1}$)', fontsize=14)
        ax.set_ylabel('Frequency (cm$^{-1}$)', fontsize=14)
        ax.set_title('Phonon Density of States', fontsize=16, pad=15)

        # 图例
        if projected_dos:
            ax.legend(loc='best', frameon=False)

        ax.grid(True, alpha=0.3, linestyle=':')

        # 保存
        if output_path:
            self.save_figure(fig, output_path)
            plt.close()
            return output_path

        plt.show()
        return None


class InteractivePhononPlotter:
    """交互式声子谱绘图器（使用Plotly）"""

    @staticmethod
    def plot_band_interactive(qpoints: np.ndarray, frequencies: np.ndarray,
                              output_path: Optional[Path] = None) -> Path:
        """
        绘制交互式声子能带图

        Parameters
        ----------
        qpoints : np.ndarray
            Q点数组
        frequencies : np.ndarray
            频率数组
        output_path : Path, optional
            输出HTML文件路径

        Returns
        -------
        Path
            输出文件路径
        """
        fig = go.Figure()

        # 添加每条能带
        nbands = frequencies.shape[1]
        for i in range(nbands):
            fig.add_trace(go.Scatter(
                x=qpoints,
                y=frequencies[:, i],
                mode='lines',
                name=f'Band {i+1}',
                showlegend=False,
                line=dict(color='blue', width=1.5)
            ))

        # 布局
        fig.update_layout(
            title='Interactive Phonon Band Structure',
            xaxis_title='Wave vector',
            yaxis_title='Frequency (cm⁻¹)',
            hovermode='closest',
            template='plotly_white',
            width=900,
            height=600
        )

        # 保存
        if output_path:
            fig.write_html(str(output_path))
            print(f"✓ 已保存交互式图表: {output_path}")
            return output_path

        fig.show()
        return None
