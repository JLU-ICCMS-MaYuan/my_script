"""
VASP Pipeline工具函数

提供Pipeline通用的辅助功能

作者：Claude
创建时间：2025-11-20
"""

import logging
from pathlib import Path
from typing import List, Dict, Any

logger = logging.getLogger(__name__)


def check_vasp_completion(work_dir: Path) -> bool:
    """
    检查VASP计算是否正常完成

    Parameters
    ----------
    work_dir : Path
        VASP计算目录

    Returns
    -------
    bool
        完成返回True
    """
    outcar = work_dir / "OUTCAR"

    if not outcar.exists():
        return False

    try:
        with open(outcar, 'r') as f:
            content = f.read()

        # 检查正常结束标志
        if "reached required accuracy" in content:
            return True

        if "writing wavefunctions" in content:
            return True

    except Exception as e:
        logger.error(f"读取OUTCAR失败: {e}")

    return False


def check_convergence(work_dir: Path, criterion: str = "energy") -> bool:
    """
    检查计算收敛性

    Parameters
    ----------
    work_dir : Path
        计算目录
    criterion : str
        收敛判据：'energy', 'force', 'stress'

    Returns
    -------
    bool
        收敛返回True
    """
    outcar = work_dir / "OUTCAR"

    if not outcar.exists():
        return False

    try:
        with open(outcar, 'r') as f:
            lines = f.readlines()

        # 查找最后的能量变化
        for line in reversed(lines):
            if "reached required accuracy" in line:
                return True

            if "WARNING" in line:
                logger.warning(f"发现WARNING: {line.strip()}")

    except Exception as e:
        logger.error(f"检查收敛性失败: {e}")

    return False


def extract_final_energy(work_dir: Path) -> float:
    """
    从OUTCAR提取最终能量

    Parameters
    ----------
    work_dir : Path
        计算目录

    Returns
    -------
    float
        能量（eV），失败返回None
    """
    outcar = work_dir / "OUTCAR"

    if not outcar.exists():
        return None

    try:
        with open(outcar, 'r') as f:
            lines = f.readlines()

        # 查找最后一个能量
        for line in reversed(lines):
            if "energy  without entropy" in line:
                parts = line.split()
                return float(parts[-1])

    except Exception as e:
        logger.error(f"提取能量失败: {e}")

    return None


def generate_summary_report(
    pipelines_results: List[Dict[str, Any]],
    output_file: Path
):
    """
    生成批量计算汇总报告

    Parameters
    ----------
    pipelines_results : List[Dict]
        所有pipeline的结果列表
    output_file : Path
        报告输出文件
    """
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("VASP Pipeline批量计算汇总报告\n")
        f.write("="*80 + "\n\n")

        total = len(pipelines_results)
        success = sum(1 for r in pipelines_results if r.get('success'))
        failed = total - success

        f.write(f"总计算数: {total}\n")
        f.write(f"成功: {success}\n")
        f.write(f"失败: {failed}\n\n")

        f.write("-"*80 + "\n")
        f.write("详细结果:\n")
        f.write("-"*80 + "\n\n")

        for i, result in enumerate(pipelines_results, 1):
            f.write(f"{i}. {result['structure_name']}\n")
            f.write(f"   工作目录: {result['work_dir']}\n")
            f.write(f"   状态: {'成功' if result.get('success') else '失败'}\n")

            if result.get('energy'):
                f.write(f"   最终能量: {result['energy']:.6f} eV\n")

            if result.get('error'):
                f.write(f"   错误信息: {result['error']}\n")

            f.write("\n")

    logger.info(f"汇总报告已生成: {output_file}")


def validate_structure_file(structure_file: Path) -> bool:
    """
    验证结构文件是否有效

    Parameters
    ----------
    structure_file : Path
        结构文件路径

    Returns
    -------
    bool
        有效返回True
    """
    if not structure_file.exists():
        logger.error(f"结构文件不存在: {structure_file}")
        return False

    if structure_file.suffix not in ['.vasp', '.POSCAR', '']:
        logger.warning(f"文件后缀不是.vasp或POSCAR: {structure_file}")

    # 简单检查文件内容
    try:
        with open(structure_file, 'r') as f:
            lines = f.readlines()

        if len(lines) < 8:
            logger.error(f"POSCAR文件行数不足: {structure_file}")
            return False

        # 检查scaling factor（第二行应该是数字）
        try:
            float(lines[1].strip())
        except:
            logger.error(f"POSCAR格式错误（第二行不是数字）: {structure_file}")
            return False

        return True

    except Exception as e:
        logger.error(f"验证结构文件失败: {e}")
        return False
