#!/usr/bin/env python
"""
VASP 命令行接口

提供现代化的命令行工具，用于执行VASP计算流程。

作者：Claude
创建时间：2025-11-20
"""

import argparse
import logging
import sys
import re
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

# 导入Pipeline和Workflow类
from vasp.pipelines import PropertiesPipeline, PhononPropertiesPipeline, BatchPipeline
from vasp.pipelines.relax import RelaxPipeline
from vasp.pipelines.md import MdPipeline
from vasp.pipelines.utils import validate_structure_file

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 默认支持的结构后缀
SUPPORTED_STRUCTURE_EXTS = ["vasp", "poscar", "cif", "res", "xsf"]


def parse_structure_exts(ext_str: Optional[str]) -> list[str]:
    if not ext_str:
        return SUPPORTED_STRUCTURE_EXTS
    parts = [p.strip().lower() for p in re.split(r"[ ,]+", ext_str) if p.strip()]
    return parts or SUPPORTED_STRUCTURE_EXTS


def parse_pressures(values) -> list[float]:
    if values is None:
        return [0.0]
    if isinstance(values, (int, float)):
        return [float(values)]
    return [float(v) for v in values]


def format_pressure_dir(pressure: float) -> str:
    return f"{pressure:g}_GPa"


def derive_work_root(input_path: Path) -> Path:
    """
    根据输入路径自动推导工作根目录：
    - 文件：使用文件所在目录下的 stem 作为工作根（去掉后缀）
    - 目录：使用目录本身作为工作根
    """
    if input_path.is_dir():
        return Path.cwd()
    return Path.cwd() / input_path.stem


MODULES = [
    "relax",
    "phonon",
    "md",
    "scf",
    "dos",
    "band",
    "elf",
    "cohp",
    "bader",
    "fermisurface",
]
PROPERTY_MODULES = {"scf", "dos", "band", "elf", "cohp", "bader", "fermisurface"}


def _rewrite_modules(argv: List[str]) -> List[str]:
    """支持多模块串写：如 `vasp relax phonon dos ...` 重写为 combo 子命令。"""
    modules: List[str] = []
    for token in argv[1:]:
        if token.startswith("-"):
            break
        if token in MODULES:
            modules.append(token)
        else:
            break
    if len(modules) > 1:
        return [argv[0], "combo", *modules, *argv[1 + len(modules):]]
    return argv


sys.argv = _rewrite_modules(sys.argv)


def run_properties_command(args, title: str, modules: list[str]):
    """通用电子性质子命令入口（scf/dos/band/elf/cohp/bader/fermisurface）。"""
    logger.info("=" * 80)
    logger.info(f"VASP {title}")
    logger.info("=" * 80)

    config: Dict[str, Any] = {}
    if args.json:
        config = load_json_config(Path(args.json))

    final_config = merge_configs(config, args)

    input_path = Path(final_config["input"])
    structure_exts = parse_structure_exts(final_config.get("structure_ext"))
    pressures = parse_pressures(final_config.get("pressure"))
    base_root = derive_work_root(input_path)
    is_batch = detect_batch_mode(input_path, structure_exts)
    tasks = final_config.get("tasks")
    parallel_flag = tasks is not None and tasks > 1
    max_workers = tasks or 1
    include_elf = "elf" in modules
    include_cohp = "cohp" in modules
    include_bader = "bader" in modules
    include_fermi = "fermisurface" in modules

    try:
        for p in pressures:
            pressure_label = format_pressure_dir(p)
            pressure_dir = base_root / pressure_label

            pipeline_kwargs = {
                "kspacing": final_config.get("kspacing", 0.2),
                "encut": final_config.get("encut"),
                "include_elf": include_elf,
                "include_cohp": include_cohp,
                "include_bader": include_bader,
                "include_fermi": include_fermi,
                "plot_dos_type": final_config.get("dos_type", "element"),
                "queue_system": final_config.get("job_system", "bash"),
                "mpi_procs": final_config.get("mpi_procs"),
                "potcar_dir": Path(final_config["potcar_dir"]) if final_config.get("potcar_dir") else None,
                "potcar_type": final_config.get("potcar_type", "PBE"),
                "prepare_only": not final_config.get("submit", False),
                "requested_steps": modules,
                "run_relax": True,
                "pressure": p,
            }

            if is_batch:
                logger.info(f"批量模式: 处理目录 {input_path}，压强 {pressures}")
                structures = scan_structure_files(input_path, structure_exts)
                if not structures:
                    logger.error("批量模式未找到结构文件")
                    sys.exit(1)

                tasks_limit = max_workers if max_workers > 0 else 1
                results = _run_matrix_tasks(
                    structures=structures,
                    pressures=pressures,
                    work_root=base_root,
                    pipeline_class=PropertiesPipeline,
                    pipeline_kwargs=pipeline_kwargs,
                    tasks_limit=tasks_limit,
                    pressure_first=False,
                )
                success_count = sum(1 for r in results if r.get("success"))
                logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")
                break
            else:
                logger.info(f"单文件模式: {input_path}, 压强 {p} GPa")
                pressure_dir.mkdir(parents=True, exist_ok=True)
                pipeline = PropertiesPipeline(
                    structure_file=input_path,
                    work_dir=pressure_dir,
                    **pipeline_kwargs,
                )
                success = pipeline.run()
                if success:
                    logger.info("\n✓ 计算完成")
                    logger.info(f"结果保存在: {pressure_dir}")
                else:
                    logger.error("\n✗ 计算失败")
                    sys.exit(1)
    except Exception as exc:
        logger.error(f"\n计算异常: {exc}", exc_info=True)
        sys.exit(1)


def load_json_config(json_file: Path) -> Dict[str, Any]:
    """
    从JSON文件加载配置

    Parameters
    ----------
    json_file : Path
        JSON配置文件路径

    Returns
    -------
    Dict[str, Any]
        配置字典
    """
    try:
        with open(json_file, 'r') as f:
            config = json.load(f)
        logger.info(f"从JSON文件加载配置: {json_file}")
        return config
    except Exception as e:
        logger.error(f"加载JSON配置文件失败: {e}")
        sys.exit(1)


def merge_configs(json_config: Dict[str, Any], cli_args: argparse.Namespace) -> Dict[str, Any]:
    """
    合并JSON配置和命令行参数

    优先级：命令行参数 > JSON配置 > 默认值

    Parameters
    ----------
    json_config : Dict[str, Any]
        JSON配置字典
    cli_args : argparse.Namespace
        命令行参数

    Returns
    -------
    Dict[str, Any]
        合并后的配置
    """
    merged = {}

    # 1. 先应用JSON配置
    for key, value in json_config.items():
        merged[key] = value

    # 2. 命令行参数覆盖JSON配置
    # 只有当命令行参数不是None（即用户明确指定）时才覆盖
    for key, value in vars(cli_args).items():
        # 跳过None值（用户未指定）和一些特殊字段
        if value is not None and key not in ['command', 'json', 'func']:
            merged[key] = value

    return merged


def detect_batch_mode(input_path: Path, structure_exts: Optional[list[str]] = None) -> bool:
    """
    智能检测是否为批量模式

    Parameters
    ----------
    input_path : Path
        输入路径
    structure_exts : list[str], optional
         结构文件后缀筛选

    Returns
    -------
    bool
        True表示批量模式
    """
    if not input_path.is_dir():
        return False

    # 检测目录中是否有结构文件
    patterns = []
    for ext in structure_exts or SUPPORTED_STRUCTURE_EXTS:
        e = ext.lower()
        if e in ("vasp", "poscar"):
            patterns.extend(['*.vasp', '*.POSCAR', 'POSCAR*'])
        elif e == "cif":
            patterns.append("*.cif")
        elif e == "res":
            patterns.append("*.res")
        elif e == "xsf":
            patterns.append("*.xsf")
        else:
            patterns.append(f"*.{e}")
    for pattern in patterns:
        files = list(input_path.glob(pattern))
        if len(files) >= 1:
            logger.info(f"检测到批量模式: 找到 {len(files)} 个结构文件")
            return True

    return False


def scan_structure_files(structures_dir: Path, structure_exts: Optional[list[str]] = None) -> list[Path]:
    """扫描目录中的结构文件。"""
    patterns = []
    for ext in structure_exts or SUPPORTED_STRUCTURE_EXTS:
        e = ext.lower()
        if e in ("vasp", "poscar"):
            patterns.extend(['*.vasp', '*.POSCAR', 'POSCAR*'])
        elif e == "cif":
            patterns.append("*.cif")
        elif e == "res":
            patterns.append("*.res")
        elif e == "xsf":
            patterns.append("*.xsf")
        else:
            patterns.append(f"*.{e}")

    files: list[Path] = []
    seen = set()
    for pattern in patterns:
        for f in structures_dir.glob(pattern):
            if f in seen:
                continue
            if validate_structure_file(f):
                files.append(f)
                seen.add(f)
    return sorted(files)


def _run_single_pipeline(
    pipeline_cls,
    structure_file: Path,
    work_dir: Path,
    pipeline_kwargs: Dict[str, Any],
    pressure: float,
):
    """运行单个 pipeline，返回结果字典。"""
    try:
        work_dir.mkdir(parents=True, exist_ok=True)
        pipeline = pipeline_cls(
            structure_file=structure_file,
            work_dir=work_dir,
            pressure=pressure,
            **pipeline_kwargs,
        )
        success = pipeline.run()
        return {"structure": structure_file, "pressure": pressure, "success": success}
    except Exception as exc:
        logger.error(f"{pipeline_cls.__name__} 执行异常: {exc}", exc_info=True)
        return {"structure": structure_file, "pressure": pressure, "success": False, "error": str(exc)}


def _run_combo_pipeline(
    modules: List[str],
    structure_file: Path,
    work_dir: Path,
    pressure: float,
    config: Dict[str, Any],
):
    """针对 combo 模式的单任务执行。"""
    try:
        work_dir.mkdir(parents=True, exist_ok=True)
        cfg_common = {
            "kspacing": config.get("kspacing", 0.2),
            "encut": config.get("encut"),
            "potcar_dir": Path(config["potcar_dir"]) if config.get("potcar_dir") else None,
            "potcar_type": config.get("potcar_type", "PBE"),
            "job_system": config.get("job_system", "bash"),
            "mpi_procs": config.get("mpi_procs"),
            "prepare_only": not config.get("submit", False),
        }

        need_phonon = "phonon" in modules
        property_modules = [m for m in modules if m in PROPERTY_MODULES]
        need_md = "md" in modules

        relax_res = _run_relax_pipeline(structure_file, work_dir, cfg_common, pressure)
        if not relax_res.get("success"):
            logger.error("relax 未完成或报错，已停止后续步骤")
            return {"structure": structure_file, "pressure": pressure, "success": False, "error": "relax failed"}

        if cfg_common.get("prepare_only", True):
            logger.info("prepare_only=True，仅生成 relax 输入，未准备后续 scf/phonon/md。")
            return {"structure": structure_file, "pressure": pressure, "success": True}

        relaxed_poscar = relax_res.get("primitive") or relax_res.get("relaxed") or structure_file
        source_structure = relaxed_poscar if Path(relaxed_poscar).exists() else structure_file

        futures = []
        success = True
        with ThreadPoolExecutor(max_workers=3) as executor:
            if need_phonon:
                futures.append(executor.submit(_run_phonon_pipeline, source_structure, work_dir, cfg_common, pressure))
            if property_modules:
                futures.append(executor.submit(_run_properties_pipeline, source_structure, work_dir, property_modules, cfg_common, pressure))
            if need_md:
                futures.append(executor.submit(_run_md_pipeline, source_structure, work_dir, cfg_common, pressure))

            for f in as_completed(futures):
                if not f.result():
                    success = False

        if not success:
            return {"structure": structure_file, "pressure": pressure, "success": False, "error": "sub pipeline failed"}
        return {"structure": structure_file, "pressure": pressure, "success": True}
    except Exception as exc:
        logger.error(f"combo 执行异常: {exc}", exc_info=True)
        return {"structure": structure_file, "pressure": pressure, "success": False, "error": str(exc)}


def _run_matrix_tasks(
    structures: List[Path],
    pressures: List[float],
    work_root: Path,
    pipeline_class,
    pipeline_kwargs: Dict[str, Any],
    tasks_limit: int,
    pressure_first: bool = False,
    combo_mode: bool = False,
):
    """
    结构×压强任务矩阵调度，tasks_limit 控制并发。
    """
    tasks: List[Tuple[Path, float]] = []
    for s in structures:
        for p in pressures:
            tasks.append((s, p))

    results = []
    max_workers = max(tasks_limit, 1)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {}
        for structure_file, pressure in tasks:
            work_dir = work_root / structure_file.stem / format_pressure_dir(pressure)
            if combo_mode:
                future = executor.submit(
                    _run_combo_pipeline,
                    pipeline_kwargs["modules"],
                    structure_file,
                    work_dir,
                    pressure,
                    pipeline_kwargs["config"],
                )
            else:
                future = executor.submit(
                    _run_single_pipeline,
                    pipeline_class,
                    structure_file,
                    work_dir,
                    pipeline_kwargs,
                    pressure,
                )
            future_map[future] = (structure_file, pressure)

        for fut in as_completed(future_map):
            res = fut.result()
            results.append(res)
            logger.info(f"完成: {res.get('structure').name} @ {res.get('pressure')} GPa -> {'成功' if res.get('success') else '失败'}")
    return results


def command_relax(args):
    """执行结构优化命令"""
    logger.info("=" * 80)
    logger.info("VASP 结构优化")
    logger.info("=" * 80)

    config = {}
    if args.json:
        config = load_json_config(Path(args.json))

    final_config = merge_configs(config, args)

    input_path = Path(final_config["input"])
    structure_exts = parse_structure_exts(final_config.get("structure_ext"))
    pressures = parse_pressures(final_config.get("pressure"))
    base_root = derive_work_root(input_path)
    is_batch = detect_batch_mode(input_path, structure_exts)
    tasks = final_config.get("tasks")
    parallel_flag = tasks is not None and tasks > 1
    max_workers = tasks or 1

    try:
        for p in pressures:
            pressure_dirname = format_pressure_dir(p)
            pipeline_kwargs = {
                "kspacing": final_config.get("kspacing", 0.2),
                "encut": final_config.get("encut"),
                "potcar_dir": Path(final_config["potcar_dir"]) if final_config.get("potcar_dir") else None,
                "potcar_type": final_config.get("potcar_type", "PBE"),
                "queue_system": final_config.get("job_system", "bash"),
                "mpi_procs": final_config.get("mpi_procs"),
                "prepare_only": not final_config.get("submit", False),
                "pressure": p,
            }

            if is_batch:
                logger.info(f"批量模式: 处理目录 {input_path}，压强 {pressures}")
                structures = scan_structure_files(input_path, structure_exts)
                if not structures:
                    logger.error("批量模式未找到结构文件")
                    sys.exit(1)

                tasks_limit = max_workers if max_workers > 0 else 1
                results = _run_matrix_tasks(
                    structures=structures,
                    pressures=pressures,
                    work_root=base_root,
                    pipeline_class=RelaxPipeline,
                    pipeline_kwargs=pipeline_kwargs,
                    tasks_limit=tasks_limit,
                    pressure_first=False,
                )
                success_count = sum(1 for r in results if r.get("success"))
                logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")
            else:
                logger.info(f"单文件模式: {input_path}, 压强 {p} GPa")
                pressure_dir = base_root / pressure_dirname
                pressure_dir.mkdir(parents=True, exist_ok=True)
                pipeline = RelaxPipeline(
                    structure_file=input_path,
                    work_dir=pressure_dir,
                    **pipeline_kwargs,
                )
                success = pipeline.run()
                if success:
                    logger.info("\n✓ 结构优化完成")
                else:
                    logger.error("\n✗ 结构优化失败")
                    sys.exit(1)
    except Exception as exc:
        logger.error(f"\n计算异常: {exc}", exc_info=True)
        sys.exit(1)


def command_scf(args):
    run_properties_command(args, "自洽计算", ["scf"])


def command_dos(args):
    run_properties_command(args, "态密度", ["dos"])


def command_band(args):
    run_properties_command(args, "能带", ["band"])


def command_elf(args):
    run_properties_command(args, "ELF", ["relax", "scf", "elf"])


def command_cohp(args):
    run_properties_command(args, "COHP", ["relax", "scf", "cohp"])


def command_bader(args):
    run_properties_command(args, "Bader 电荷", ["relax", "scf", "bader"])


def command_fermi(args):
    run_properties_command(args, "费米面", ["relax", "scf", "fermisurface"])


def command_properties(args):
    """统一入口，便于子命令共享参数。"""
    modules = getattr(args, "modules", None) or getattr(args, "base_steps", ["scf"])
    title = getattr(args, "title", "电子性质")
    run_properties_command(args, title, modules)


def _attach_common_property_args(parser: argparse.ArgumentParser):
    parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    parser.add_argument('--json', help='JSON配置文件路径')
    parser.add_argument('--tasks', type=int, help='同时运行的最大结构数（并行度，默认串行）')
    parser.add_argument('--kspacing', type=float, help='K点间距')
    parser.add_argument('--encut', type=float, help='截断能(eV)')
    parser.add_argument('--dos-type', choices=['element', 'spd', 'element_spd'], help='DOS投影类型')
    parser.add_argument('--potcar-dir', help='POTCAR库目录')
    parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    parser.add_argument('-p', '--pressure', type=float, nargs='+', help='外压(GPa)，可多值')
    parser.add_argument('--structure-ext', type=str, help='目录输入时的结构后缀过滤，逗号分隔，默认vasp')
    parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    parser.add_argument('--mpi-procs', type=str, help="MPI 启动命令，可为数字(默认 mpirun -np N)或完整前缀，如 'mpirun -np 16' / 'srun -n 16'")
    parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')

def command_phonon(args):
    """执行声子性质计算命令"""
    logger.info("=" * 80)
    logger.info("VASP 声子性质全流程计算")
    logger.info("=" * 80)

    # 加载JSON配置（如果有）
    config = {}
    if args.json:
        config = load_json_config(Path(args.json))

    # 合并配置
    final_config = merge_configs(config, args)

    input_path = Path(final_config['input'])
    structure_exts = parse_structure_exts(final_config.get('structure_ext'))
    pressures = parse_pressures(final_config.get('pressure'))
    base_root = derive_work_root(input_path)
    # 检测批量模式
    is_batch = detect_batch_mode(input_path, structure_exts)
    tasks = final_config.get("tasks")
    parallel_flag = tasks is not None and tasks > 1
    max_workers = tasks or 1

    try:
        pipeline_kwargs = {
            'supercell': final_config.get('supercell', [2, 2, 2]),
            'method': final_config.get('method', 'disp'),
            'kspacing': final_config.get('kspacing', 0.3),
            'encut': final_config.get('encut'),
            'queue_system': final_config.get('job_system', 'bash'),
            'mpi_procs': final_config.get('mpi_procs'),
            'potcar_dir': Path(final_config['potcar_dir']) if final_config.get('potcar_dir') else None,
            'potcar_type': final_config.get('potcar_type', 'PBE'),
            'prepare_only': not final_config.get('submit', False),
            'include_relax': True,
        }

        if is_batch:
            logger.info(f"批量模式: 处理目录 {input_path}，压强 {pressures}")
            structures = scan_structure_files(input_path, structure_exts)
            if not structures:
                logger.error("批量模式未找到结构文件")
                sys.exit(1)

            tasks_limit = max_workers if max_workers > 0 else 1
            results = _run_matrix_tasks(
                structures=structures,
                pressures=pressures,
                work_root=base_root,
                pipeline_class=PhononPropertiesPipeline,
                pipeline_kwargs=pipeline_kwargs,
                tasks_limit=tasks_limit,
                pressure_first=False,
            )
            success_count = sum(1 for r in results if r.get('success'))
            logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")

        else:
            for p in pressures:
                pressure_label = format_pressure_dir(p)
                pressure_dir = base_root / pressure_label
                pressure_dir.mkdir(parents=True, exist_ok=True)

                pipeline = PhononPropertiesPipeline(
                    structure_file=input_path,
                    work_dir=pressure_dir,
                    pressure=p,
                    **pipeline_kwargs
                )

                success = pipeline.run()

                if success:
                    logger.info("\n✓ 计算成功完成！")
                    logger.info(f"结果保存在: {pressure_dir}")
                else:
                    logger.error("\n✗ 计算失败")
                    sys.exit(1)

    except Exception as e:
        logger.error(f"\n计算异常: {e}", exc_info=True)
        sys.exit(1)


def command_md(args):
    """执行分子动力学命令"""
    logger.info("=" * 80)
    logger.info("VASP 分子动力学")
    logger.info("=" * 80)

    config = {}
    if args.json:
        config = load_json_config(Path(args.json))

    final_config = merge_configs(config, args)

    input_path = Path(final_config["input"])
    structure_exts = parse_structure_exts(final_config.get("structure_ext"))
    pressures = parse_pressures(final_config.get("pressure"))
    base_root = derive_work_root(input_path)
    is_batch = detect_batch_mode(input_path, structure_exts)
    tasks = final_config.get("tasks")
    parallel_flag = tasks is not None and tasks > 1
    max_workers = tasks or 1

    try:
        for p in pressures:
            pressure_label = format_pressure_dir(p)
            pressure_dir = base_root / pressure_label

            pipeline_kwargs = {
                "potim": final_config.get("potim", 1.0),
                "tebeg": final_config.get("tebeg", 300.0),
                "teend": final_config.get("teend", 300.0),
                "nsw": final_config.get("nsw", 200),
                "kspacing": final_config.get("kspacing", 0.2),
                "encut": final_config.get("encut"),
                "potcar_dir": Path(final_config["potcar_dir"]) if final_config.get("potcar_dir") else None,
                "potcar_type": final_config.get("potcar_type", "PBE"),
                "queue_system": final_config.get("job_system", "bash"),
                "mpi_procs": final_config.get("mpi_procs"),
                "prepare_only": not final_config.get("submit", False),
                "include_relax": True,
                "pressure": p,
            }

            if is_batch:
                logger.info(f"批量模式: 压强 {p} GPa, 处理目录 {input_path}")
                batch = BatchPipeline(
                    pipeline_class=MdPipeline,
                    structures_dir=input_path,
                    work_root=base_root,
                    pipeline_kwargs=pipeline_kwargs,
                    parallel=parallel_flag,
                    max_workers=max_workers,
                    structure_exts=structure_exts,
                    pressure_label=pressure_label,
                )
                results = batch.run()
                success_count = sum(1 for r in results if r.get('success'))
                logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功 (压强 {p} GPa)")
            else:
                pressure_dir.mkdir(parents=True, exist_ok=True)
                pipeline = MdPipeline(
                    structure_file=input_path,
                    work_dir=pressure_dir,
                    **pipeline_kwargs,
                )
                success = pipeline.run()

                if success:
                    logger.info("\n✓ 分子动力学计算完成")
                    logger.info(f"结果保存在: {pressure_dir}")
                else:
                    logger.error("\n✗ 分子动力学计算失败")
                    sys.exit(1)
    except Exception as exc:
        logger.error(f"\n计算异常: {exc}", exc_info=True)
        sys.exit(1)


def _run_relax_pipeline(structure_file: Path, work_dir: Path, cfg: Dict[str, Any], pressure: float):
    pipeline = RelaxPipeline(
        structure_file=structure_file,
        work_dir=work_dir,
        kspacing=cfg.get("kspacing", 0.2),
        encut=cfg.get("encut"),
        potcar_dir=cfg.get("potcar_dir"),
        potcar_type=cfg.get("potcar_type", "PBE"),
        queue_system=cfg.get("job_system", "bash"),
        mpi_procs=cfg.get("mpi_procs"),
        prepare_only=cfg.get("prepare_only", True),
        pressure=pressure,
    )
    ok = pipeline.run()
    return {
        "success": ok,
        "relaxed": work_dir / "POSCAR_relaxed",
        "primitive": Path(pipeline.steps_data.get("primitive_structure")) if pipeline.steps_data.get("primitive_structure") else None,
        "conventional": Path(pipeline.steps_data.get("conventional_structure")) if pipeline.steps_data.get("conventional_structure") else None,
    }


def _run_phonon_pipeline(structure_file: Path, work_dir: Path, cfg: Dict[str, Any], pressure: float):
    pipeline = PhononPropertiesPipeline(
        structure_file=structure_file,
        work_dir=work_dir,
        supercell=cfg.get("supercell", [2, 2, 2]),
        method=cfg.get("method", "disp"),
        kspacing=cfg.get("kspacing", 0.3),
        encut=cfg.get("encut"),
        queue_system=cfg.get("job_system", "bash"),
        mpi_procs=cfg.get("mpi_procs"),
        potcar_dir=cfg.get("potcar_dir"),
        potcar_type=cfg.get("potcar_type", "PBE"),
        prepare_only=cfg.get("prepare_only", True),
        include_relax=False,
        pressure=pressure,
        checkpoint_file=work_dir / "phonon_checkpoint.json",
        report_file=work_dir / "phonon_report.txt",
    )
    return pipeline.run()


def _run_properties_pipeline(structure_file: Path, work_dir: Path, modules: list[str], cfg: Dict[str, Any], pressure: float):
    include_elf = "elf" in modules
    include_cohp = "cohp" in modules
    include_bader = "bader" in modules
    include_fermi = "fermisurface" in modules
    pipeline = PropertiesPipeline(
        structure_file=structure_file,
        work_dir=work_dir,
        kspacing=cfg.get("kspacing", 0.2),
        encut=cfg.get("encut"),
        include_elf=include_elf,
        include_cohp=include_cohp,
        include_bader=include_bader,
        include_fermi=include_fermi,
        plot_dos_type=cfg.get("dos_type", "element"),
        queue_system=cfg.get("job_system", "bash"),
        mpi_procs=cfg.get("mpi_procs"),
        potcar_dir=cfg.get("potcar_dir"),
        potcar_type=cfg.get("potcar_type", "PBE"),
        prepare_only=cfg.get("prepare_only", True),
        requested_steps=modules,
        run_relax=False,
        pressure=pressure,
        checkpoint_file=work_dir / "electronic_checkpoint.json",
        report_file=work_dir / "electronic_report.txt",
    )
    return pipeline.run()


def _run_md_pipeline(structure_file: Path, work_dir: Path, cfg: Dict[str, Any], pressure: float):
    pipeline = MdPipeline(
        structure_file=structure_file,
        work_dir=work_dir,
        potim=cfg.get("potim", 1.0),
        tebeg=cfg.get("tebeg", 300.0),
        teend=cfg.get("teend", 300.0),
        nsw=cfg.get("nsw", 200),
        kspacing=cfg.get("kspacing", 0.2),
        encut=cfg.get("encut"),
        potcar_dir=cfg.get("potcar_dir"),
        potcar_type=cfg.get("potcar_type", "PBE"),
        queue_system=cfg.get("job_system", "bash"),
        mpi_procs=cfg.get("mpi_procs"),
        prepare_only=cfg.get("prepare_only", True),
        include_relax=False,
        pressure=pressure,
        checkpoint_file=work_dir / "md_checkpoint.json",
        report_file=work_dir / "md_report.txt",
    )
    return pipeline.run()


def command_combo(args):
    """多模块组合命令：如 `vasp relax phonon dos ...`。"""
    modules = args.modules
    logger.info("=" * 80)
    logger.info(f"VASP 组合计算: {' '.join(modules)}")
    logger.info("=" * 80)

    config = {}
    if args.json:
        config = load_json_config(Path(args.json))
    final_config = merge_configs(config, args)

    input_path = Path(final_config["input"])
    structure_exts = parse_structure_exts(final_config.get("structure_ext"))
    pressures = parse_pressures(final_config.get("pressure"))
    base_root = derive_work_root(input_path)
    is_batch = detect_batch_mode(input_path, structure_exts)
    tasks = final_config.get("tasks") or 1
    tasks = tasks if tasks > 0 else 1

    if not modules:
        logger.error("未指定任何模块")
        sys.exit(1)

    property_modules = [m for m in modules if m in PROPERTY_MODULES]
    need_phonon = "phonon" in modules
    need_md = "md" in modules
    need_relax = True  # 所有依赖 relax

    # 结构列表
    structures: list[Path]
    if is_batch:
        structures = scan_structure_files(input_path, structure_exts)
        if not structures:
            logger.error("批量模式未找到结构文件")
            sys.exit(1)
    else:
        structures = [input_path]

    try:
        logger.info(f"批量模式: 处理目录 {input_path}，压强 {pressures}")
        tasks_limit = tasks if tasks > 0 else 1
        results = _run_matrix_tasks(
            structures=structures,
            pressures=pressures,
            work_root=base_root,
            pipeline_class=None,  # combo 特殊处理
            pipeline_kwargs={
                "modules": modules,
                "need_phonon": need_phonon,
                "property_modules": property_modules,
                "need_md": need_md,
                "config": final_config,
            },
            tasks_limit=tasks_limit,
            pressure_first=False,
            combo_mode=True,
        )
        success_count = sum(1 for r in results if r.get("success"))
        logger.info(f"\n批量组合计算完成: {success_count}/{len(results)} 成功")

    except Exception as exc:
        logger.error(f"\n计算异常: {exc}", exc_info=True)
        sys.exit(1)


def create_parser():
    """创建命令行解析器"""
    parser = argparse.ArgumentParser(
        prog='vasp',
        description='VASP计算命令行工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 结构优化
  vasp relax -i POSCAR -j bash --mpi-procs 8 --submit

  # 自洽/态密度/能带
  vasp scf  -i POSCAR --submit
  vasp dos  -i POSCAR --submit
  vasp band -i POSCAR --submit

  # ELF/COHP/Bader/费米面
  vasp elf           -i POSCAR --submit
  vasp cohp          -i POSCAR --submit
  vasp bader         -i POSCAR --submit
  vasp fermisurface  -i POSCAR --submit

  # 声子/MD（默认先做 relax）
  vasp phonon -i POSCAR --supercell 2 2 2 --submit
  vasp md     -i POSCAR --potim 1 --nsw 200 --submit
  vasp combo relax phonon dos -i POSCAR --submit
        """
    )

    subparsers = parser.add_subparsers(dest='command', required=True, help='子命令')

    # ========== combo 子命令 ==========
    combo_parser = subparsers.add_parser('combo', help='多模块组合执行（可串写 relax/phonon/scf/dos/...）')
    combo_parser.add_argument('modules', nargs='+', choices=MODULES, help='模块列表，按书写顺序执行，自动补齐依赖并并行调度')
    combo_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    combo_parser.add_argument('--json', help='JSON配置文件路径')
    combo_parser.add_argument('--tasks', type=int, help='同时运行的最大结构数（默认串行）')
    combo_parser.add_argument('--kspacing', type=float, help='K点间距')
    combo_parser.add_argument('--encut', type=float, help='截断能(eV)')
    combo_parser.add_argument('--supercell', nargs=3, type=int, metavar=('X', 'Y', 'Z'), help='声子超胞，如 2 2 2')
    combo_parser.add_argument('--method', choices=['disp', 'dfpt'], help='声子计算方法')
    combo_parser.add_argument('--potim', type=float, help='MD 时间步(fs)')
    combo_parser.add_argument('--tebeg', type=float, help='MD 起始温度(K)')
    combo_parser.add_argument('--teend', type=float, help='MD 结束温度(K)')
    combo_parser.add_argument('--nsw', type=int, help='MD 步数')
    combo_parser.add_argument('--dos-type', choices=['element', 'spd', 'element_spd'], help='DOS 投影类型')
    combo_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    combo_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    combo_parser.add_argument('-p', '--pressure', type=float, nargs='+', help='外压(GPa)，可多值')
    combo_parser.add_argument('--structure-ext', type=str, help='目录输入时的结构后缀过滤，逗号分隔，默认vasp')
    combo_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    combo_parser.add_argument('--mpi-procs', type=str, help="MPI 启动命令，可为数字(默认 mpirun -np N)或完整前缀，如 'mpirun -np 16' / 'srun -n 16'")
    combo_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    combo_parser.set_defaults(func=command_combo)

    # ========== relax 子命令 ==========
    relax_parser = subparsers.add_parser('relax', help='结构优化')
    relax_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    relax_parser.add_argument('--json', help='JSON配置文件路径')
    relax_parser.add_argument('--tasks', type=int, help='同时运行的最大结构数（并行度，默认串行）')
    relax_parser.add_argument('--kspacing', type=float, help='K点间距')
    relax_parser.add_argument('--encut', type=float, help='截断能(eV)')
    relax_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    relax_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    relax_parser.add_argument('-p', '--pressure', type=float, nargs='+', help='外压(GPa)，可多值')
    relax_parser.add_argument('--structure-ext', type=str, help='目录输入时的结构后缀过滤，逗号分隔，默认vasp')
    relax_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    relax_parser.add_argument('--mpi-procs', type=str, help="MPI 启动命令，可为数字(默认 mpirun -np N)或完整前缀，如 'mpirun -np 16' / 'srun -n 16'")
    relax_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    relax_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    relax_parser.set_defaults(func=command_relax)

    # ========== 性质类子命令（scf/dos/band/elf/cohp/bader/fermisurface） ==========
    property_commands = [
        ("scf", "自洽计算", ["scf"]),
        ("dos", "态密度", ["dos"]),
        ("band", "能带", ["band"]),
        ("elf", "ELF", ["elf"]),
        ("cohp", "COHP", ["cohp"]),
        ("bader", "Bader 电荷", ["bader"]),
        ("fermisurface", "费米面", ["fermisurface"]),
    ]
    for name, title, modules in property_commands:
        prop_parser = subparsers.add_parser(name, help=title)
        _attach_common_property_args(prop_parser)
        prop_parser.set_defaults(func=command_properties, modules=modules, title=title)

    # ========== phonon 子命令 ==========
    phonon_parser = subparsers.add_parser('phonon', help='声子性质全流程计算')
    phonon_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    phonon_parser.add_argument('--json', help='JSON配置文件路径')
    phonon_parser.add_argument('--tasks', type=int, help='同时运行的最大结构数（并行度，默认串行）')
    phonon_parser.add_argument('--supercell', nargs=3, type=int, metavar=('X', 'Y', 'Z'), help='超胞大小，如: --supercell 2 2 2')
    phonon_parser.add_argument('--method', choices=['disp', 'dfpt'], help='声子计算方法')
    phonon_parser.add_argument('--kspacing', type=float, help='K点间距')
    phonon_parser.add_argument('--encut', type=float, help='截断能(eV)')
    phonon_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    phonon_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    phonon_parser.add_argument('-p', '--pressure', type=float, nargs='+', help='外压(GPa)，可多值')
    phonon_parser.add_argument('--structure-ext', type=str, help='目录输入时的结构后缀过滤，逗号分隔，默认vasp')
    phonon_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    phonon_parser.add_argument('--mpi-procs', type=str, help="MPI 启动命令，可为数字(默认 mpirun -np N)或完整前缀，如 'mpirun -np 16' / 'srun -n 16'")
    phonon_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    phonon_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    phonon_parser.set_defaults(func=command_phonon)

    # ========== md 子命令 ==========
    md_parser = subparsers.add_parser('md', help='分子动力学')
    md_parser.add_argument('-i', '--input', required=True, help='输入文件')
    md_parser.add_argument('--json', help='JSON配置文件路径')
    md_parser.add_argument('--potim', type=float, help='时间步长(fs)')
    md_parser.add_argument('--tebeg', type=float, help='起始温度(K)')
    md_parser.add_argument('--teend', type=float, help='结束温度(K)')
    md_parser.add_argument('--nsw', type=int, help='MD步数')
    md_parser.add_argument('--kspacing', type=float, help='K点间距')
    md_parser.add_argument('--encut', type=float, help='截断能(eV)')
    md_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    md_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    md_parser.add_argument('-p', '--pressure', type=float, nargs='+', help='外压(GPa)，可多值')
    md_parser.add_argument('--structure-ext', type=str, help='目录输入时的结构后缀过滤，逗号分隔，默认vasp')
    md_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    md_parser.add_argument('--mpi-procs', type=str, help="MPI 启动命令，可为数字(默认 mpirun -np N)或完整前缀，如 'mpirun -np 16' / 'srun -n 16'")
    md_parser.add_argument('--tasks', type=int, help='同时运行的最大结构数（并行度，默认串行）')
    md_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    md_parser.set_defaults(func=command_md)

    return parser


def main():
    """主入口函数"""
    parser = create_parser()
    args = parser.parse_args()

    # 设置日志级别
    if hasattr(args, 'log_level') and args.log_level:
        logging.getLogger().setLevel(getattr(logging, args.log_level))

    # 执行子命令
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
