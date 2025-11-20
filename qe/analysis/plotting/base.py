"""
绘图基类模块

提供发表级质量的绘图功能，支持matplotlib和plotly。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Optional, Dict, Any, List
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import matplotlib as mpl


# 设置发表级质量的matplotlib参数
def setup_publication_style():
    """配置发表级质量的绘图样式"""
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['figure.dpi'] = 300  # 高DPI
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['savefig.bbox'] = 'tight'
    mpl.rcParams['legend.frameon'] = False


class BasePlotter(ABC):
    """绘图基类"""

    def __init__(self, style: str = 'publication'):
        """
        初始化绘图器

        Parameters
        ----------
        style : str
            绘图风格 ('publication', 'presentation', 'simple')
        """
        self.style = style
        if style == 'publication':
            setup_publication_style()

    @abstractmethod
    def plot(self, data: Any, output_path: Optional[Path] = None, **kwargs) -> Path:
        """
        绘图方法（子类实现）

        Parameters
        ----------
        data : Any
            绘图数据
        output_path : Path, optional
            输出文件路径
        **kwargs
            其他参数

        Returns
        -------
        Path
            输出文件路径
        """
        pass

    def save_figure(self, fig: plt.Figure, output_path: Path, formats: List[str] = ['png', 'pdf']):
        """
        保存图片（多种格式）

        Parameters
        ----------
        fig : Figure
            matplotlib图形对象
        output_path : Path
            输出路径（不含扩展名）
        formats : List[str]
            输出格式列表
        """
        saved_files = []
        for fmt in formats:
            file_path = output_path.with_suffix(f'.{fmt}')
            fig.savefig(file_path, format=fmt, dpi=300, bbox_inches='tight')
            saved_files.append(file_path)
            print(f"✓ 已保存: {file_path}")

        return saved_files
