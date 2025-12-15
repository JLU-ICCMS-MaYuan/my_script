"""
VASP Pipeline工具函数

提供Pipeline通用的辅助功能

作者：Claude
创建时间：2025-11-20
"""

import logging
import os
import shutil
import math
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from ase.io import read as ase_read, write as ase_write
import spglib

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
    # 仅检查存在，其余解析交给 ensure_poscar/ASE 处理
    return True


def ensure_poscar(src: Path, dest: Path) -> Path:
    """
    将任意受支持的结构文件转换/复制为 POSCAR。

    Parameters
    ----------
    src : Path
        原始结构文件路径
    dest : Path
        目标 POSCAR 路径

    Returns
    -------
    Path
        生成的 POSCAR 路径
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        # POSCAR/vasp 直接复制
        if src.suffix.lower() in [".vasp", ".poscar"] or src.name.upper().startswith("POSCAR"):
            shutil.copy(src, dest)
        else:
            atoms = ase_read(src)
            ase_write(dest, atoms, format="vasp")
    except Exception as exc:
        logger.error(f"转换结构为 POSCAR 失败: {src} -> {dest}，错误: {exc}")
        raise
    return dest


def find_symmetry(poscar: Path, output_dir: Path, symprec: float = 1e-3) -> tuple[Optional[Path], Optional[Path], Optional[str]]:
    """
    使用 spglib 生成原胞和标准晶胞，并返回路径与空间群。
    """
    try:
        atoms = ase_read(poscar)
        cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        spacegroup = spglib.get_spacegroup(cell, symprec=symprec)

        prim_lat, prim_pos, prim_num = spglib.find_primitive(cell, symprec=symprec)
        std_lat, std_pos, std_num = spglib.standardize_cell(cell, symprec=symprec)

        prim_atoms = ase_read(poscar)
        prim_atoms.set_cell(prim_lat, scale_atoms=False)
        prim_atoms.set_scaled_positions(prim_pos)

        std_atoms = ase_read(poscar)
        std_atoms.set_cell(std_lat, scale_atoms=False)
        std_atoms.set_scaled_positions(std_pos)

        output_dir.mkdir(parents=True, exist_ok=True)
        prim_path = output_dir / "POSCAR_primitive"
        std_path = output_dir / "POSCAR_conventional"
        ase_write(prim_path, prim_atoms, format="vasp", direct=True)
        ase_write(std_path, std_atoms, format="vasp", direct=True)

        info_file = output_dir / "symmetry.txt"
        info_file.write_text(f"Spacegroup: {spacegroup}\nSymprec: {symprec}\nPrimitive: {prim_path.name}\nConventional: {std_path.name}\n")

        return prim_path, std_path, spacegroup
    except Exception as exc:
        logger.warning(f"对称性分析失败: {exc}")
        return None, None, None


def prepare_potcar(
    poscar_file: Path,
    potcar_dir: Optional[Path],
    output_file: Path,
    potcar_type: str = "PBE",
    potcar_lib_root: Optional[Path] = None,
) -> bool:
    """
    准备POTCAR文件（用户优先选择）。
    优先使用当前工作目录下的 potcar_lib 中已选赝势；若缺失则从 potcar_dir 搜索唯一匹配的元素目录复制到 potcar_lib。
    若找到多个候选或找不到候选，提示用户手动放置到 potcar_lib 后重试。

    Parameters
    ----------
    poscar_file : Path
        POSCAR文件路径
    potcar_dir : Path
        POTCAR源库目录（包含各元素的POTCAR文件），可为 None（仅使用现有 potcar_lib）
    output_file : Path
        输出的POTCAR文件路径
    potcar_type : str
        POTCAR类型：'PBE', 'LDA', 'PW91'等
    potcar_lib_root : Path, optional
        赝势缓存目录，默认使用当前工作目录下的 potcar_lib

    Returns
    -------
    bool
        成功返回True
    """
    try:
        elements = _parse_poscar_elements(poscar_file)
        if not elements:
            logger.error("未能从 POSCAR 解析元素列表")
            return False

        logger.info(f"从POSCAR读取到元素: {elements}")

        potcar_lib = potcar_lib_root or Path.cwd() / "potcar_lib"
        potcar_lib.mkdir(parents=True, exist_ok=True)

        potcar_content: List[str] = []

        for element in elements:
            # 1) 先看 potcar_lib 是否已有选定的赝势
            existing = _locate_potcar_in_lib(potcar_lib, element)
            if existing:
                logger.info(f"使用 potcar_lib 中的赝势: {existing}")
                potcar_content.append(existing.read_text())
                continue

            # 2) 从 potcar_dir 搜索候选
            if not potcar_dir:
                logger.error(f"potcar_lib 缺少 {element}，且未提供 potcar_dir，请将所选 POTCAR 放入 {potcar_lib} 后重试")
                return False

            candidates = _search_potcar_candidates(potcar_dir, element, potcar_type)
            if len(candidates) == 0:
                logger.error(f"未找到元素 {element} 的POTCAR，请在 {potcar_lib} 手动放置或检查 potcar_dir")
                return False

            chosen = _prompt_choose_potcar(element, candidates)
            if not chosen:
                logger.error(f"未为元素 {element} 选择赝势，已中止")
                return False

            dest = potcar_lib / element
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(chosen, dest)
            logger.info(f"已将 {element} 赝势复制到 {dest}，后续重复计算将复用")
            potcar_content.append(Path(dest).read_text())

        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text("".join(potcar_content))
        logger.info(f"POTCAR已生成: {output_file}")
        return True

    except Exception as e:
        logger.error(f"准备POTCAR失败: {e}", exc_info=True)
        return False


def _parse_poscar_elements(poscar_file: Path) -> List[str]:
    """简单解析 POSCAR 的元素行（第6行），返回元素符号列表。"""
    try:
        lines = poscar_file.read_text().splitlines()
        if len(lines) < 6:
            return []
        candidates = lines[5].split()
        # 如果这一行是数字，说明不含符号，直接失败
        if all(item.replace('.', '', 1).isdigit() for item in candidates):
            return []
        return candidates
    except Exception:
        return []


def _locate_potcar_in_lib(potcar_lib: Path, element: str) -> Optional[Path]:
    """
    在 potcar_lib 中查找指定元素的 POTCAR。
    支持两种布局：
      - potcar_lib/Element  (文件)
      - potcar_lib/Element/POTCAR  (目录+文件)
    """
    file_candidate = potcar_lib / element
    dir_candidate = potcar_lib / element / "POTCAR"
    if dir_candidate.exists():
        return dir_candidate
    if file_candidate.exists():
        return file_candidate
    return None


def _search_potcar_candidates(potcar_dir: Path, element: str, potcar_type: str) -> List[Path]:
    """
    在 potcar_dir 下搜索元素的 POTCAR 候选。
    匹配规则：目录名以元素符号开头，内部有 POTCAR。
    """
    candidates: List[Path] = []
    base_dir = Path(potcar_dir)

    # 允许 potcar_type/Element/POTCAR 或 Element*/POTCAR
    search_roots = [base_dir, base_dir / potcar_type] if potcar_type else [base_dir]
    for root in search_roots:
        if not root.exists():
            continue
        for child in root.iterdir():
            if child.is_dir() and child.name.lower().startswith(element.lower()):
                pot_path = child / "POTCAR"
                if pot_path.exists():
                    candidates.append(pot_path)
    return candidates


def _prompt_choose_potcar(element: str, candidates: List[Path]) -> Optional[Path]:
    """
    交互式选择赝势候选。
    """
    if len(candidates) == 1:
        return candidates[0]

    print(f"\n为元素 {element} 选择赝势（按序号输入，回车默认1）：")
    for idx, c in enumerate(candidates, 1):
        print(f"  {idx}. {c}")
    choice = input("请输入编号：").strip()
    if choice == "":
        return candidates[0]
    try:
        num = int(choice)
        if 1 <= num <= len(candidates):
            return candidates[num - 1]
    except Exception:
        pass
    logger.error("输入无效，未选择赝势")
    return None


def kspacing_to_mesh(poscar_file: Path, kspacing: float) -> Tuple[int, int, int]:
    """
    将 KSPACING 转换为 Monkhorst-Pack 网格：N_i = max(1, ceil(|b_i| / kspacing))
    其中 b_i 为 2π 倒格矢长度。
    """
    try:
        lines = poscar_file.read_text().splitlines()
        if len(lines) < 5:
            raise ValueError("POSCAR 行数不足")
        scale = float(lines[1].split()[0])
        a1 = [float(x) for x in lines[2].split()[:3]]
        a2 = [float(x) for x in lines[3].split()[:3]]
        a3 = [float(x) for x in lines[4].split()[:3]]
        a1 = [v * scale for v in a1]
        a2 = [v * scale for v in a2]
        a3 = [v * scale for v in a3]

        def cross(u, v):
            return [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]

        def dot(u, v):
            return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

        vol = dot(a1, cross(a2, a3))
        if abs(vol) < 1e-8:
            raise ValueError("晶胞体积异常")

        factor = 2 * math.pi / vol
        b1 = [x * factor for x in cross(a2, a3)]
        b2 = [x * factor for x in cross(a3, a1)]
        b3 = [x * factor for x in cross(a1, a2)]

        def norm(u):
            return math.sqrt(dot(u, u))

        n1 = max(1, math.ceil(norm(b1) / kspacing))
        n2 = max(1, math.ceil(norm(b2) / kspacing))
        n3 = max(1, math.ceil(norm(b3) / kspacing))
        return int(n1), int(n2), int(n3)
    except Exception as exc:
        logger.error(f"KSPACING 转换 KPOINTS 失败: {exc}")
        # 回退默认 1x1x1，避免写入错误浮点
        return 1, 1, 1
