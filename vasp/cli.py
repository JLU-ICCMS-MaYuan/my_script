#!/usr/bin/env python
"""
VASP 命令行接口

提供现代化的命令行工具，用于执行VASP计算流程。

用法示例：
    # 单个结构的电子性质计算
    vasp electronic -i POSCAR -w ./calc --kspacing 0.2 --include-elf

    # 批量计算（自动检测）
    vasp electronic -i ./structures/ -w ./batch_calc --kspacing 0.2

    # 使用JSON配置文件
    vasp electronic -i POSCAR -w ./calc --json config.json

    # 声子计算
    vasp phonon -i POSCAR -w ./calc --supercell 2 2 2

作者：Claude
创建时间：2025-11-20
"""

import argparse
import json
import logging
import sys
import re
from pathlib import Path
from typing import Dict, Any, Optional

# 导入Pipeline和Workflow类
from vasp.pipelines import ElectronicPropertiesPipeline, PhononPropertiesPipeline, BatchPipeline
from vasp.pipelines.relax import RelaxPipeline
from vasp.pipelines.md import MdPipeline
from vasp.utils.job import submit_job

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def _rewrite_relax_combo(argv: list[str]) -> list[str]:
    """
    支持“vasp relax phonon ...”或“vasp relax md ...”语法：
    等价于“vasp phonon --with-relax ...”或“vasp md --with-relax ...”。
    仅当前两个非选项参数分别为 relax 和 {phonon,md} 时生效。
    """
    if len(argv) >= 3:
        first, second = argv[1], argv[2]
        if first == "relax" and second in {"phonon", "md"}:
            return [argv[0], second, "--with-relax", *argv[3:]]
    return argv


sys.argv = _rewrite_relax_combo(sys.argv)


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


def detect_batch_mode(input_path: Path, batch_flag: bool) -> bool:
    """
    智能检测是否为批量模式

    Parameters
    ----------
    input_path : Path
        输入路径
    batch_flag : bool
        用户是否明确指定--batch

    Returns
    -------
    bool
        True表示批量模式
    """
    # 用户明确指定
    if batch_flag:
        return True

    # 输入不是目录
    if not input_path.is_dir():
        return False

    # 检测目录中是否有多个结构文件
    patterns = ['*.vasp', '*.POSCAR', 'POSCAR*']
    for pattern in patterns:
        files = list(input_path.glob(pattern))
        if len(files) >= 1:  # 只要有结构文件就认为是批量
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
    work_dir = Path(final_config["work_dir"])
    is_batch = detect_batch_mode(input_path, final_config.get("batch", False))

    pipeline_kwargs = {
        "kspacing": final_config.get("kspacing", 0.2),
        "encut": final_config.get("encut"),
        "potcar_dir": Path(final_config["potcar_dir"]) if final_config.get("potcar_dir") else None,
        "potcar_type": final_config.get("potcar_type", "PBE"),
        "queue_system": final_config.get("job_system", "bash"),
        "mpi_procs": final_config.get("mpi_procs"),
        "submit_only": final_config.get("submit", False),
        "prepare_only": not final_config.get("submit", False),
    }

    try:
        if is_batch:
            logger.info(f"批量模式: 处理目录 {input_path}")

            batch = BatchPipeline(
                pipeline_class=RelaxPipeline,
                structures_dir=input_path,
                work_root=work_dir,
                pipeline_kwargs=pipeline_kwargs,
                parallel=final_config.get("parallel", False),
                max_workers=final_config.get("max_workers", 4),
            )
            results = batch.run()
            success_count = sum(1 for r in results if r.get("success"))
            logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")
        else:
            logger.info(f"单文件模式: {input_path}")
            pipeline = RelaxPipeline(
                structure_file=input_path,
                work_dir=work_dir,
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


def command_electronic(args):
    """执行电子性质计算命令"""
    logger.info("=" * 80)
    logger.info("VASP 电子性质全流程计算")
    logger.info("=" * 80)

    # 加载JSON配置（如果有）
    config = {}
    if args.json:
        config = load_json_config(Path(args.json))

    # 合并配置
    final_config = merge_configs(config, args)

    input_path = Path(final_config['input'])
    work_dir = Path(final_config['work_dir'])

    # 检测批量模式
    is_batch = detect_batch_mode(input_path, final_config.get('batch', False))

    master_script = final_config.get('master_script', False)
    submit_flag = final_config.get('submit', False)
    steps_raw = final_config.get('steps')
    custom_steps = None
    if steps_raw:
        parts = re.split(r'[ ,]+', steps_raw)
        custom_steps = [p for p in parts if p]

    # 准备Pipeline参数
    pipeline_kwargs = {
        'kspacing': final_config.get('kspacing', 0.2),
        'encut': final_config.get('encut'),
        'include_elf': final_config.get('include_elf', False),
        'include_cohp': final_config.get('include_cohp', False),
        'plot_dos_type': final_config.get('dos_type', 'element'),
        'queue_system': final_config.get('job_system', 'bash'),
        'mpi_procs': final_config.get('mpi_procs'),
        'potcar_dir': Path(final_config['potcar_dir']) if final_config.get('potcar_dir') else None,
        'potcar_type': final_config.get('potcar_type', 'PBE'),
        'submit_only': submit_flag if not master_script else False,
        'prepare_only': (not submit_flag) or master_script,
        'custom_steps': custom_steps,
        'master_script': master_script,
        'master_script_name': final_config.get('master_script_name', 'run_master.sh'),
    }

    try:
        if is_batch:
            # 批量计算模式
            logger.info(f"批量模式: 处理目录 {input_path}")

            batch = BatchPipeline(
                pipeline_class=ElectronicPropertiesPipeline,
                structures_dir=input_path,
                work_root=work_dir,
                pipeline_kwargs=pipeline_kwargs,
                parallel=final_config.get('parallel', False),
                max_workers=final_config.get('max_workers', 4),
            )

            results = batch.run()

            # 打印结果摘要
            success_count = sum(1 for r in results if r.get('success'))
            logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")
            if master_script and submit_flag:
                logger.info("批量 + --master-script 当前仅生成各目录的 run_master.sh，不会自动提交，请按需手动提交。")

        else:
            # 单个文件模式
            logger.info(f"单文件模式: {input_path}")

            pipeline = ElectronicPropertiesPipeline(
                structure_file=input_path,
                work_dir=work_dir,
                **pipeline_kwargs
            )

            success = pipeline.run()

            if success:
                if master_script:
                    script_path = pipeline.master_script_path or (work_dir / final_config.get('master_script_name', 'run_master.sh'))
                    if submit_flag:
                        if not script_path.exists():
                            logger.error(f"未找到总控脚本: {script_path}")
                            sys.exit(1)
                        job_id = submit_job(script_path, final_config.get('job_system', 'bash'))
                        logger.info(f"已提交总控脚本: {job_id}")
                    else:
                        logger.info(f"已生成总控脚本（未提交）: {script_path}")
                logger.info("\n✓ 计算成功完成！")
                logger.info(f"结果保存在: {work_dir}")
            else:
                logger.error("\n✗ 计算失败")
                sys.exit(1)

    except Exception as e:
        logger.error(f"\n计算异常: {e}", exc_info=True)
        sys.exit(1)


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
    work_dir = Path(final_config['work_dir'])

    # 检测批量模式
    is_batch = detect_batch_mode(input_path, final_config.get('batch', False))

    steps_raw = final_config.get('steps')
    custom_steps = None
    if steps_raw:
        custom_steps = [p for p in re.split(r'[ ,]+', steps_raw) if p]

    # 准备Pipeline参数
    pipeline_kwargs = {
        'supercell': final_config.get('supercell', [2, 2, 2]),
        'method': final_config.get('method', 'disp'),
        'kspacing': final_config.get('kspacing', 0.3),
        'encut': final_config.get('encut'),
        'queue_system': final_config.get('job_system', 'bash'),
        'mpi_procs': final_config.get('mpi_procs'),
        'potcar_dir': Path(final_config['potcar_dir']) if final_config.get('potcar_dir') else None,
        'potcar_type': final_config.get('potcar_type', 'PBE'),
        'submit_only': final_config.get('submit', False),
        'prepare_only': not final_config.get('submit', False),
        'include_relax': final_config.get('with_relax', False),
        'custom_steps': custom_steps,
    }

    try:
        if is_batch:
            # 批量计算模式
            logger.info(f"批量模式: 处理目录 {input_path}")

            batch = BatchPipeline(
                pipeline_class=PhononPropertiesPipeline,
                structures_dir=input_path,
                work_root=work_dir,
                pipeline_kwargs=pipeline_kwargs,
                parallel=final_config.get('parallel', False),
                max_workers=final_config.get('max_workers', 4),
            )

            results = batch.run()

            success_count = sum(1 for r in results if r.get('success'))
            logger.info(f"\n批量计算完成: {success_count}/{len(results)} 成功")

        else:
            # 单个文件模式
            logger.info(f"单文件模式: {input_path}")

            pipeline = PhononPropertiesPipeline(
                structure_file=input_path,
                work_dir=work_dir,
                **pipeline_kwargs
            )

            success = pipeline.run()

            if success:
                logger.info("\n✓ 计算成功完成！")
                logger.info(f"结果保存在: {work_dir}")
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
    work_dir = Path(final_config["work_dir"])

    steps_raw = final_config.get("steps")
    custom_steps = [p for p in re.split(r"[ ,]+", steps_raw) if p] if steps_raw else None

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
        "submit_only": final_config.get("submit", False),
        "prepare_only": not final_config.get("submit", False),
        "include_relax": final_config.get("with_relax", False),
        "custom_steps": custom_steps,
    }

    try:
        pipeline = MdPipeline(
            structure_file=input_path,
            work_dir=work_dir,
            **pipeline_kwargs,
        )
        success = pipeline.run()

        if success:
            logger.info("\n✓ 分子动力学计算完成")
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
  # 电子性质计算
  vasp electronic -i POSCAR -w ./calc --kspacing 0.2 --include-elf

  # 批量计算
  vasp electronic -i ./structures/ -w ./batch_calc --batch --parallel

  # 使用JSON配置
  vasp electronic -i POSCAR -w ./calc --json config.json

  # 声子计算
  vasp phonon -i POSCAR -w ./calc --supercell 2 2 2
        """
    )

    subparsers = parser.add_subparsers(dest='command', required=True, help='子命令')

    # ========== relax 子命令 ==========
    relax_parser = subparsers.add_parser('relax', help='结构优化')
    relax_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    relax_parser.add_argument('-w', '--work-dir', required=True, help='工作目录')
    relax_parser.add_argument('--json', help='JSON配置文件路径')
    relax_parser.add_argument('--batch', action='store_true', help='批量模式')
    relax_parser.add_argument('--parallel', action='store_true', help='并行执行批量任务')
    relax_parser.add_argument('--max-workers', type=int, help='最大并行数（默认4）')
    relax_parser.add_argument('--kspacing', type=float, help='K点间距')
    relax_parser.add_argument('--encut', type=float, help='截断能(eV)')
    relax_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    relax_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    relax_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    relax_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    relax_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    relax_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    relax_parser.set_defaults(func=command_relax)

    # ========== electronic 子命令 ==========
    electronic_parser = subparsers.add_parser('electronic', help='电子性质全流程计算')
    electronic_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    electronic_parser.add_argument('-w', '--work-dir', required=True, help='工作目录')
    electronic_parser.add_argument('--json', help='JSON配置文件路径')
    electronic_parser.add_argument('--batch', action='store_true', help='批量模式')
    electronic_parser.add_argument('--parallel', action='store_true', help='并行执行批量任务')
    electronic_parser.add_argument('--max-workers', type=int, help='最大并行数')
    electronic_parser.add_argument('--kspacing', type=float, help='K点间距')
    electronic_parser.add_argument('--encut', type=float, help='截断能(eV)')
    electronic_parser.add_argument('--include-elf', action='store_true', help='包含ELF计算')
    electronic_parser.add_argument('--include-cohp', action='store_true', help='包含COHP计算')
    electronic_parser.add_argument('--dos-type', choices=['element', 'spd', 'element_spd'], help='DOS投影类型')
    electronic_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    electronic_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    electronic_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    electronic_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    electronic_parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,scf,dos,band,elf,cohp,plotting')
    electronic_parser.add_argument('--master-script', action='store_true', help='生成总控脚本，一次提交串行跑完所选步骤')
    electronic_parser.add_argument('--master-script-name', type=str, default='run_master.sh', help='总控脚本文件名')
    electronic_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    electronic_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    electronic_parser.set_defaults(func=command_electronic)

    # ========== phonon 子命令 ==========
    phonon_parser = subparsers.add_parser('phonon', help='声子性质全流程计算')
    phonon_parser.add_argument('-i', '--input', required=True, help='输入文件或目录')
    phonon_parser.add_argument('-w', '--work-dir', required=True, help='工作目录')
    phonon_parser.add_argument('--json', help='JSON配置文件路径')
    phonon_parser.add_argument('--batch', action='store_true', help='批量模式')
    phonon_parser.add_argument('--parallel', action='store_true', help='并行执行批量任务')
    phonon_parser.add_argument('--max-workers', type=int, help='最大并行数')
    phonon_parser.add_argument('--supercell', nargs=3, type=int, metavar=('X', 'Y', 'Z'), help='超胞大小，如: --supercell 2 2 2')
    phonon_parser.add_argument('--method', choices=['disp', 'dfpt'], help='声子计算方法')
    phonon_parser.add_argument('--kspacing', type=float, help='K点间距')
    phonon_parser.add_argument('--encut', type=float, help='截断能(eV)')
    phonon_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    phonon_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    phonon_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    phonon_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    phonon_parser.add_argument('--with-relax', action='store_true', help='先做结构优化再做声子')
    phonon_parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,phonon_prepare,phonon_calculate,phonon_band,phonon_dos,plotting')
    phonon_parser.add_argument('--submit', action='store_true', help='提交作业（默认仅生成输入和脚本）')
    phonon_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'], help='日志级别')
    phonon_parser.set_defaults(func=command_phonon)

    # ========== md 子命令 ==========
    md_parser = subparsers.add_parser('md', help='分子动力学')
    md_parser.add_argument('-i', '--input', required=True, help='输入文件')
    md_parser.add_argument('-w', '--work-dir', required=True, help='工作目录')
    md_parser.add_argument('--json', help='JSON配置文件路径')
    md_parser.add_argument('--potim', type=float, help='时间步长(fs)')
    md_parser.add_argument('--tebeg', type=float, help='起始温度(K)')
    md_parser.add_argument('--teend', type=float, help='结束温度(K)')
    md_parser.add_argument('--nsw', type=int, help='MD步数')
    md_parser.add_argument('--kspacing', type=float, help='K点间距')
    md_parser.add_argument('--encut', type=float, help='截断能(eV)')
    md_parser.add_argument('--potcar-dir', help='POTCAR库目录')
    md_parser.add_argument('--potcar-type', choices=['PBE', 'LDA', 'PW91'], help='POTCAR类型')
    md_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs', 'lsf'], help='队列系统')
    md_parser.add_argument('--mpi-procs', type=int, help='MPI进程数（默认取配置或8）')
    md_parser.add_argument('--with-relax', action='store_true', help='先做结构优化再进行MD')
    md_parser.add_argument('--steps', type=str, help='自定义步骤序列，逗号分隔: relax,md')
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
