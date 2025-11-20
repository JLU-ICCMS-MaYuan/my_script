#!/usr/bin/env python
"""
EPW 命令行接口

提供现代化的命令行工具，用于执行EPW计算流程。

用法示例：
    # EPW完整计算流程
    epw run -i relax.out -w ./calc --mode full

    # 仅电声耦合计算
    epw run -i relax.out -w ./calc --mode elph --nkf 20 20 20

    # 使用JSON配置
    epw run -i relax.out -w ./calc --json config.json

作者：Claude
创建时间：2025-11-20
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Any, Optional

# 导入Workflow类
# TODO: 等workflows/__init__.py修复后再启用
# from epw.workflows import EPWWorkflow

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_json_config(json_file: Path) -> Dict[str, Any]:
    """从JSON文件加载配置"""
    try:
        with open(json_file, 'r') as f:
            config = json.load(f)
        logger.info(f"从JSON文件加载配置: {json_file}")
        return config
    except Exception as e:
        logger.error(f"加载JSON配置文件失败: {e}")
        sys.exit(1)


def merge_configs(json_config: Dict[str, Any], cli_args: argparse.Namespace) -> Dict[str, Any]:
    """合并JSON配置和命令行参数"""
    merged = {}

    # 1. 先应用JSON配置
    for key, value in json_config.items():
        merged[key] = value

    # 2. 命令行参数覆盖JSON配置
    for key, value in vars(cli_args).items():
        if value is not None and key not in ['command', 'json', 'func']:
            merged[key] = value

    return merged


def command_run(args):
    """执行EPW计算命令"""
    logger.info("=" * 80)
    logger.info("EPW 电声耦合计算")
    logger.info("=" * 80)

    # 加载JSON配置（如果有）
    config = {}
    if args.json:
        config = load_json_config(Path(args.json))

    # 合并配置
    final_config = merge_configs(config, args)

    # TODO: 实现EPWWorkflow调用
    logger.info("EPW计算功能开发中...")
    logger.info(f"输入: {final_config['input']}")
    logger.info(f"工作目录: {final_config['work_dir']}")
    logger.info(f"模式: {final_config.get('mode', 'full')}")
    logger.info(f"配置: {final_config}")


def create_parser():
    """创建命令行解析器"""
    parser = argparse.ArgumentParser(
        prog='epw',
        description='EPW电声耦合计算命令行工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # EPW完整计算流程
  epw run -i relax.out -w ./calc --mode full

  # 仅声子插值
  epw run -i relax.out -w ./calc --mode phono --nqf 20 20 20

  # 超导性质计算
  epw run -i relax.out -w ./calc --mode sc

  # 使用JSON配置
  epw run -i relax.out -w ./calc --json config.json
        """
    )

    subparsers = parser.add_subparsers(dest='command', required=True, help='子命令')

    # ========== run 子命令 ==========
    run_parser = subparsers.add_parser('run', help='EPW计算')
    run_parser.add_argument('-i', '--input', required=True, help='输入文件(relax.out等)')
    run_parser.add_argument('-w', '--work-dir', required=True, help='工作目录')
    run_parser.add_argument('--json', help='JSON配置文件路径')
    run_parser.add_argument('--mode',
                           choices=['full', 'eband', 'phono', 'elph', 'sc'],
                           help='计算模式')
    run_parser.add_argument('--nkf', nargs=3, type=int, metavar=('X', 'Y', 'Z'),
                           help='精细k点网格')
    run_parser.add_argument('--nqf', nargs=3, type=int, metavar=('X', 'Y', 'Z'),
                           help='精细q点网格')
    run_parser.add_argument('--pp-dir', help='QE赝势库目录')
    run_parser.add_argument('-j', '--job-system', choices=['bash', 'slurm', 'pbs'],
                           help='队列系统')
    run_parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING'],
                           help='日志级别')
    run_parser.set_defaults(func=command_run)

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
