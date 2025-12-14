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
from typing import Dict, Any, Optional
import json

# 导入Pipeline和Workflow类
from vasp.pipelines import PropertiesPipeline, PhononPropertiesPipeline, BatchPipeline
from vasp.pipelines.relax import RelaxPipeline
from vasp.pipelines.md import MdPipeline

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_structure_exts(ext_str: Optional[str]) -> list[str]:
    if not ext_str:
        return ["vasp"]
    parts = [p.strip().lower() for p in re.split(r"[ ,]+", ext_str) if p.strip()]
    return parts or ["vasp"]


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


PROPERTY_ALLOWED_STEPS = {
    "relax",
    "scf",
    "dos",
    "band",
    "elf",
    "cohp",
    "bader",
    "fermisurface",
    "plotting",
}


def _rewrite_relax_combo(argv: list[str]) -> list[str]:
    """
    支持“vasp relax phonon ...”或“vasp relax md ...”等语法：
    仅当前两个非选项参数分别为 relax 和受支持的子命令时，将其重写为对应子命令。
    """
    if len(argv) >= 3:
        first, second = argv[1], argv[2]
        supported = {
            "phonon",
            "md",
            "scf",
            "dos",
            "band",
            "elf",
            "cohp",
            "bader",
            "fermisurface",
        }
        if first == "relax" and second in supported:
            return [argv[0], second, *argv[3:]]
    return argv


sys.argv = _rewrite_relax_combo(sys.argv)


def normalize_property_steps(base_steps: list[str], override: Optional[str]) -> list[str]:
    """
    合并默认步骤与用户自定义步骤，并自动补齐依赖关系。
    """
    raw_steps = base_steps[:]
    if override:
        raw_steps = [p for p in re.split(r"[ ,]+", override) if p]

    deps = {
        "scf": ["relax"],
        "dos": ["relax", "scf"],
        "band": ["relax", "scf"],
        "elf": ["relax", "scf"],
        "cohp": ["relax", "scf"],
        "bader": ["relax", "scf"],
        "fermisurface": ["relax", "scf"],
    }

    ordered: list[str] = []

    def add_step(step: str):
        name = step.strip().lower()
        if not name:
            return
        if name not in PROPERTY_ALLOWED_STEPS:
            logger.warning(f"忽略未支持的步骤: {name}")
            return
        for dep in deps.get(name, []):
            add_step(dep)
        if name not in ordered:
            ordered.append(name)

    for step in raw_steps:
        add_step(step)

    return ordered


def run_properties_command(args, title: str, base_steps: list[str]):
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
    steps = normalize_property_steps(base_steps, final_config.get("steps"))

    if not steps:
        logger.error("未获得有效步骤序列，请检查 --steps 参数")
        sys.exit(1)

    include_elf = "elf" in steps
    include_cohp = "cohp" in steps
    include_bader = "bader" in steps
    include_fermi = "fermisurface" in steps

    try:
        for p in pressures:
            pressure_dir = base_root / format_pressure_dir(p)
            pressure_dir.mkdir(parents=True, exist_ok=True)

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
                "custom_steps": steps,
                "pressure": p,
            }

            if is_batch:
                logger.info(f"批量模式: 压强 {p} GPa, 处理目录 {input_path}")
                batch = BatchPipeline(
                    pipeline_class=PropertiesPipeline,
                    structures_dir=input_path,
                    work_root=base_root,
                    pipeline_kwargs=pipeline_kwargs,
                    parallel=parallel_flag,
                    max_workers=max_workers,
                    structure_exts=structure_exts,
                    pressure_label=format_pressure_dir(p),
                )
                results = batch.run()
                success_count = sum(1 for r in results if r.get("success"))
                logger.info(f"\n批量计算完成(压强 {p} GPa): {success_count}/{len(results)} 成功")
            else:
                logger.info(f"单文件模式: {input_path}, 压强 {p} GPa")
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
    for ext in structure_exts or ["vasp"]:
        e = ext.lower()
        if e == "vasp":
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
                logger.info(f"批量模式: 压强 {p} GPa, 处理目录 {input_path}")

                batch = BatchPipeline(
                    pipeline_class=RelaxPipeline,
                    structures_dir=input_path,
                    work_root=base_root,
                    pipeline_kwargs=pipeline_kwargs,
                    parallel=parallel_flag,
                    max_workers=max_workers,
                    structure_exts=structure_exts,
                    pressure_label=pressure_dirname,
                )
                results = batch.run()
                success_count = sum(1 for r in results if r.get("success"))
                logger.info(f"\n批量计算完成(压强 {p} GPa): {success_count}/{len(results)} 成功")
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
    run_properties_command(args, "自洽计算", ["relax", "scf"])


def command_dos(args):
    run_properties_command(args, "态密度", ["relax", "scf", "dos"])


def command_band(args):
    run_properties_command(args, "能带", ["relax", "scf", "band"])


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
    base_steps = getattr(args, "base_steps", ["relax", "scf"])
    title = getattr(args, "title", "电子性质")
    run_properties_command(args, title, base_steps)


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
    parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,scf,dos,band,elf,cohp,bader,fermisurface')
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

    steps_raw = final_config.get('steps')
    custom_steps = None
    if steps_raw:
        custom_steps = [p for p in re.split(r'[ ,]+', steps_raw) if p]
        if "relax" not in [s.lower() for s in custom_steps]:
            custom_steps.insert(0, "relax")

    try:
        for p in pressures:
            pressure_dir = base_root / format_pressure_dir(p)
            pressure_dir.mkdir(parents=True, exist_ok=True)

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
                'custom_steps': custom_steps,
                'pressure': p,
            }

            if is_batch:
                # 批量计算模式
                logger.info(f"批量模式: 压强 {p} GPa, 处理目录 {input_path}")

                batch = BatchPipeline(
                    pipeline_class=PhononPropertiesPipeline,
                    structures_dir=input_path,
                    work_root=base_root,
                    pipeline_kwargs=pipeline_kwargs,
                    parallel=parallel_flag,
                    max_workers=max_workers,
                    structure_exts=structure_exts,
                    pressure_label=pressure_dir.name,
                )

                results = batch.run()

                success_count = sum(1 for r in results if r.get('success'))
                logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功 (压强 {p} GPa)")

            else:
                # 单个文件模式
                logger.info(f"单文件模式: {input_path}, 压强 {p} GPa")

                pipeline = PhononPropertiesPipeline(
                    structure_file=input_path,
                    work_dir=pressure_dir,
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

    steps_raw = final_config.get("steps")
    custom_steps = [p for p in re.split(r"[ ,]+", steps_raw) if p] if steps_raw else None
    if custom_steps and "relax" not in [s.lower() for s in custom_steps]:
        custom_steps.insert(0, "relax")

    try:
        for p in pressures:
            pressure_dir = base_root / format_pressure_dir(p)
            pressure_dir.mkdir(parents=True, exist_ok=True)

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
                "custom_steps": custom_steps,
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
                    pressure_label=pressure_dir.name,
                )
                results = batch.run()
                success_count = sum(1 for r in results if r.get('success'))
                logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功 (压强 {p} GPa)")
            else:
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
  vasp dos  -i POSCAR --steps relax,scf,dos,band --submit
  vasp band -i POSCAR --submit

  # ELF/COHP/Bader/费米面
  vasp elf           -i POSCAR --submit
  vasp cohp          -i POSCAR --submit
  vasp bader         -i POSCAR --submit
  vasp fermisurface  -i POSCAR --submit

  # 声子/MD（默认先做 relax）
  vasp phonon -i POSCAR --supercell 2 2 2 --submit
  vasp md     -i POSCAR --potim 1 --nsw 200 --submit
        """
    )

    subparsers = parser.add_subparsers(dest='command', required=True, help='子命令')

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
    relax_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    relax_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    relax_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    relax_parser.set_defaults(func=command_relax)

    # ========== 性质类子命令（scf/dos/band/elf/cohp/bader/fermisurface） ==========
    property_commands = [
        ("scf", "自洽计算", ["relax", "scf"]),
        ("dos", "态密度", ["relax", "scf", "dos"]),
        ("band", "能带", ["relax", "scf", "band"]),
        ("elf", "ELF", ["relax", "scf", "elf"]),
        ("cohp", "COHP", ["relax", "scf", "cohp"]),
        ("bader", "Bader 电荷", ["relax", "scf", "bader"]),
        ("fermisurface", "费米面", ["relax", "scf", "fermisurface"]),
    ]
    for name, title, steps in property_commands:
        prop_parser = subparsers.add_parser(name, help=title)
        _attach_common_property_args(prop_parser)
        prop_parser.set_defaults(func=command_properties, base_steps=steps, title=title)

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
    phonon_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    phonon_parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,phonon_prepare,phonon_calculate,phonon_band,phonon_dos,plotting')
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
    md_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    md_parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,md')
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
