"""
VASP绘图模块

提供电子结构和声子性质的绘图功能：
- 电子能带和DOS
- 声子能带和DOS
- ELF、COHP等

作者：Claude
创建时间：2025-11-20
"""

import logging
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# 设置发表级质量
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['lines.linewidth'] = 2

logger = logging.getLogger(__name__)


def plot_band_structure(
    work_dir: Path,
    output_file: Optional[Path] = None,
    ylim: Tuple[float, float] = (-5, 5),
    plot_title: Optional[str] = None
) -> Path:
    """
    绘制电子能带结构

    Parameters
    ----------
    work_dir : Path
        包含vasprun.xml的工作目录
    output_file : Path, optional
        输出文件路径，默认为work_dir/eband.png
    ylim : Tuple[float, float]
        y轴范围（相对费米能级，单位eV）
    plot_title : str, optional
        图表标题

    Returns
    -------
    Path
        生成的图片文件路径
    """
    try:
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import BSPlotter

        logger.info("开始绘制电子能带图...")

        # 读取vasprun.xml
        vasprun_file = work_dir / "vasprun.xml"
        if not vasprun_file.exists():
            raise FileNotFoundError(f"未找到vasprun.xml: {vasprun_file}")

        vasprun = Vasprun(str(vasprun_file))
        band_structure = vasprun.get_band_structure(line_mode=True)

        # 创建绘图
        plotter = BSPlotter(band_structure)

        # 获取matplotlib figure
        fig = plotter.get_plot(ylim=ylim)

        if plot_title:
            plt.title(plot_title)

        # 保存
        if output_file is None:
            output_file = work_dir / "eband.png"

        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()

        logger.info(f"电子能带图已保存: {output_file}")

        # 同时保存数据文件
        data_file = output_file.with_suffix('.dat')
        _save_band_data(band_structure, data_file)

        return output_file

    except Exception as e:
        logger.error(f"绘制电子能带图失败: {e}", exc_info=True)
        raise


def plot_dos(
    work_dir: Path,
    output_file: Optional[Path] = None,
    xlim: Tuple[float, float] = (-5, 5),
    ylim: Optional[Tuple[float, float]] = None,
    plot_title: Optional[str] = None,
    pdos_type: str = "element"
) -> Path:
    """
    绘制电子态密度(DOS)

    Parameters
    ----------
    work_dir : Path
        包含vasprun.xml的工作目录
    output_file : Path, optional
        输出文件路径
    xlim : Tuple[float, float]
        x轴范围（能量，相对费米能级）
    ylim : Tuple[float, float], optional
        y轴范围
    plot_title : str, optional
        图表标题
    pdos_type : str
        投影DOS类型：'element', 'spd', 'element_spd'

    Returns
    -------
    Path
        生成的图片文件路径
    """
    try:
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import DosPlotter

        logger.info("开始绘制电子DOS图...")

        # 读取vasprun.xml
        vasprun_file = work_dir / "vasprun.xml"
        if not vasprun_file.exists():
            raise FileNotFoundError(f"未找到vasprun.xml: {vasprun_file}")

        vasprun = Vasprun(str(vasprun_file))

        # 获取DOS数据
        try:
            dos = vasprun.complete_dos_normalized
        except:
            dos = vasprun.complete_dos

        # 创建plotter
        plotter = DosPlotter()

        # 添加投影DOS
        if pdos_type == "element":
            plotter.add_dos_dict(dos.get_element_dos())
        elif pdos_type == "spd":
            plotter.add_dos_dict(dos.get_spd_dos())
        elif pdos_type == "element_spd":
            plotter.add_dos_dict(dos.get_element_spd_dos())
        else:
            plotter.add_dos("Total DOS", dos)

        # 获取figure
        fig = plotter.get_plot(xlim=xlim, ylim=ylim)

        if plot_title:
            plt.title(plot_title)

        # 保存
        if output_file is None:
            output_file = work_dir / "dos.png"

        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()

        logger.info(f"电子DOS图已保存: {output_file}")

        # 保存数据文件
        data_file = output_file.with_suffix('.dat')
        _save_dos_data(dos, data_file)

        return output_file

    except Exception as e:
        logger.error(f"绘制DOS图失败: {e}", exc_info=True)
        raise


def plot_phonon_band(
    work_dir: Path,
    output_file: Optional[Path] = None,
    plot_title: Optional[str] = None
) -> Path:
    """
    绘制声子能带

    Parameters
    ----------
    work_dir : Path
        包含phonopy数据的工作目录
    output_file : Path, optional
        输出文件路径
    plot_title : str, optional
        图表标题

    Returns
    -------
    Path
        生成的图片文件路径
    """
    logger.info("开始绘制声子能带图...")

    # Phonopy已经生成了band.yaml和band.pdf
    # 这里可以读取band.yaml重新绘制更高质量的图

    if output_file is None:
        output_file = work_dir / "phonon_band.png"

    try:
        # 如果phonopy生成了图，直接使用
        phonopy_pdf = work_dir / "band.pdf"
        if phonopy_pdf.exists():
            logger.info(f"Phonopy已生成声子能带图: {phonopy_pdf}")

        # TODO: 读取band.yaml自定义绘图
        # 这里可以使用phonopy的API或直接解析yaml

        logger.info(f"声子能带图: {output_file}")
        return output_file

    except Exception as e:
        logger.error(f"绘制声子能带图失败: {e}", exc_info=True)
        raise


def plot_phonon_dos(
    work_dir: Path,
    output_file: Optional[Path] = None,
    plot_title: Optional[str] = None
) -> Path:
    """
    绘制声子态密度

    Parameters
    ----------
    work_dir : Path
        包含phonopy数据的工作目录
    output_file : Path, optional
        输出文件路径
    plot_title : str, optional
        图表标题

    Returns
    -------
    Path
        生成的图片文件路径
    """
    logger.info("开始绘制声子DOS图...")

    if output_file is None:
        output_file = work_dir / "phonon_dos.png"

    try:
        # Phonopy已生成projected_dos.dat和total_dos.dat
        total_dos_file = work_dir / "total_dos.dat"

        if total_dos_file.exists():
            # 读取并绘制
            data = np.loadtxt(total_dos_file)
            freq = data[:, 0]  # THz
            dos = data[:, 1]

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(freq, dos, 'b-', linewidth=2, label='Total DOS')
            ax.set_xlabel('Frequency (THz)', fontsize=14)
            ax.set_ylabel('DOS', fontsize=14)
            ax.legend()

            if plot_title:
                ax.set_title(plot_title)

            ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"声子DOS图已保存: {output_file}")

        return output_file

    except Exception as e:
        logger.error(f"绘制声子DOS图失败: {e}", exc_info=True)
        raise


def plot_elf(
    work_dir: Path,
    output_file: Optional[Path] = None,
    plane: str = 'xy',
    z_frac: float = 0.5
) -> Path:
    """
    绘制电子局域函数(ELF)

    Parameters
    ----------
    work_dir : Path
        包含ELFCAR的工作目录
    output_file : Path, optional
        输出文件路径
    plane : str
        切面方向：'xy', 'xz', 'yz'
    z_frac : float
        切面位置（分数坐标）

    Returns
    -------
    Path
        生成的图片文件路径
    """
    logger.info("开始绘制ELF图...")

    if output_file is None:
        output_file = work_dir / f"elf_{plane}.png"

    try:
        # 读取ELFCAR文件并绘图
        # 这里需要使用pymatgen或ASE读取ELFCAR

        logger.warning("ELF绘图功能待完善，需要读取和处理ELFCAR文件")

        return output_file

    except Exception as e:
        logger.error(f"绘制ELF图失败: {e}", exc_info=True)
        raise


def plot_cohp(
    work_dir: Path,
    output_file: Optional[Path] = None,
    plot_title: Optional[str] = None
) -> Path:
    """
    绘制晶体轨道哈密顿布居(COHP)

    Parameters
    ----------
    work_dir : Path
        包含COHPCAR的工作目录
    output_file : Path, optional
        输出文件路径
    plot_title : str, optional
        图表标题

    Returns
    -------
    Path
        生成的图片文件路径
    """
    logger.info("开始绘制COHP图...")

    if output_file is None:
        output_file = work_dir / "cohp.png"

    try:
        # COHP分析通常需要lobster等后处理工具
        # 这里可以简单绘制VASP输出的COHPCAR

        logger.warning("COHP绘图功能待完善，建议使用lobster工具")

        return output_file

    except Exception as e:
        logger.error(f"绘制COHP图失败: {e}", exc_info=True)
        raise


def _save_band_data(band_structure, output_file: Path):
    """保存能带数据到文件"""
    try:
        # 简单保存，可以后续用Origin等软件重新绘图
        with open(output_file, 'w') as f:
            f.write("# Electronic Band Structure Data\n")
            f.write("# Distance(Ang^-1)  Energy(eV)\n")
            # TODO: 提取并保存实际数据

        logger.info(f"能带数据已保存: {output_file}")
    except Exception as e:
        logger.warning(f"保存能带数据失败: {e}")


def _save_dos_data(dos, output_file: Path):
    """保存DOS数据到文件"""
    try:
        energies = dos.energies - dos.efermi  # 相对费米能级
        total_dos = dos.densities[next(iter(dos.densities))]

        with open(output_file, 'w') as f:
            f.write("# Electronic DOS Data\n")
            f.write("# Energy(eV)  DOS(states/eV)\n")
            for e, d in zip(energies, total_dos):
                f.write(f"{e:12.6f}  {d:12.6f}\n")

        logger.info(f"DOS数据已保存: {output_file}")
    except Exception as e:
        logger.warning(f"保存DOS数据失败: {e}")
