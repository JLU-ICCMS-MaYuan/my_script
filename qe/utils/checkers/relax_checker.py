"""
QE任务完成状态检查器

检查SCF、Relax等QE计算任务的完成状态。
递归搜索目录，检查输出文件中的"JOB DONE"标记。

重构自: mytoolkit/qe_toolkit/cqopt.py

作者：Claude
创建时间：2025-11-20
"""

import os
from pathlib import Path
from typing import List, Dict, Optional
from enum import Enum
from dataclasses import dataclass

from rich.console import Console
from rich.table import Table
from rich import box


class TaskStatus(str, Enum):
    """任务状态"""
    SUCCEEDED = "succeeded"  # 成功完成
    FAILED = "failed"        # 失败
    NONE = "none"            # 未找到输出文件


@dataclass
class TaskResult:
    """任务检查结果"""
    task_type: str  # relax, scf, nscf等
    path: Path
    status: TaskStatus
    output_file: str


class QETaskChecker:
    """QE任务检查器"""

    def __init__(self, base_dir: Path = None, max_depth: Optional[int] = None):
        """
        初始化检查器

        Parameters
        ----------
        base_dir : Path, optional
            基础目录，默认当前目录
        max_depth : int, optional
            最大搜索深度，None表示无限制
        """
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()
        self.max_depth = max_depth
        self.console = Console()

    def _walk_with_depth(self):
        """带深度限制的目录遍历"""
        base_dir = self.base_dir.resolve()

        for root, dirs, files in os.walk(base_dir):
            rel_path = Path(root).relative_to(base_dir)

            # 检查深度
            if self.max_depth is not None and len(rel_path.parts) > self.max_depth:
                dirs[:] = []  # 停止更深层搜索
                continue

            yield Path(root), files

    def _check_job_done(self, file_path: Path) -> bool:
        """
        检查文件中是否有JOB DONE标记

        Parameters
        ----------
        file_path : Path
            输出文件路径

        Returns
        -------
        bool
            是否完成
        """
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            return 'JOB DONE' in content
        except Exception:
            return False

    def check_tasks(self, output_filename: str) -> Dict[str, List[Path]]:
        """
        检查特定输出文件的任务状态

        Parameters
        ----------
        output_filename : str
            输出文件名（如relax.out, scf.out）

        Returns
        -------
        Dict[str, List[Path]]
            按状态分类的路径列表
        """
        results = {
            "succeeded": [],
            "failed": [],
            "none": []
        }

        for root_path, files in self._walk_with_depth():
            if output_filename in files:
                file_path = root_path / output_filename

                if self._check_job_done(file_path):
                    results["succeeded"].append(root_path)
                else:
                    results["failed"].append(root_path)
            else:
                results["none"].append(root_path)

        return results

    def check_all_types(self, output_types: List[str] = None) -> Dict[str, Dict[str, List[Path]]]:
        """
        检查多种类型的任务

        Parameters
        ----------
        output_types : List[str], optional
            要检查的输出文件类型列表，默认['relax.out', 'scf.out']

        Returns
        -------
        Dict[str, Dict[str, List[Path]]]
            {任务类型: {状态: [路径列表]}}
        """
        if output_types is None:
            output_types = ['relax.out', 'scf.out']

        all_results = {}

        for output_file in output_types:
            task_type = output_file.replace('.out', '')
            all_results[task_type] = self.check_tasks(output_file)

        return all_results

    def print_summary(self, results: Dict[str, List[Path]], task_type: str):
        """
        打印任务摘要

        Parameters
        ----------
        results : Dict[str, List[Path]]
            检查结果
        task_type : str
            任务类型
        """
        table = Table(
            title=f"{task_type.upper()} 任务状态",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold magenta",
        )

        table.add_column("状态", style="cyan", width=15)
        table.add_column("数量", style="yellow", width=10, justify="right")

        succeeded_count = len(results['succeeded'])
        failed_count = len(results['failed'])
        none_count = len(results['none'])
        total = succeeded_count + failed_count + none_count

        # 计算成功率
        if (succeeded_count + failed_count) > 0:
            success_rate = succeeded_count / (succeeded_count + failed_count) * 100
        else:
            success_rate = 0.0

        table.add_row("✓ 成功", f"[green]{succeeded_count}[/]")
        table.add_row("✗ 失败", f"[red]{failed_count}[/]")
        table.add_row("○ 无输出", f"[dim]{none_count}[/]")
        table.add_row("[bold]总计", f"[bold]{total}[/]")
        table.add_row("[bold]成功率", f"[bold green]{success_rate:.1f}%[/]")

        self.console.print(table)

    def save_results(self, results: Dict[str, List[Path]], task_type: str):
        """
        保存检查结果到文件

        Parameters
        ----------
        results : Dict[str, List[Path]]
            检查结果
        task_type : str
            任务类型
        """
        for status, paths in results.items():
            filename = f"{task_type}-{status}.txt"
            filepath = self.base_dir / filename

            with open(filepath, 'w') as f:
                for path in paths:
                    f.write(str(path.resolve()) + '\n')

            self.console.print(f"已保存: {filepath} ({len(paths)} 个路径)")

    def print_failed_details(self, results: Dict[str, List[Path]], task_type: str, max_show: int = 10):
        """
        打印失败任务的详细信息

        Parameters
        ----------
        results : Dict[str, List[Path]]
            检查结果
        task_type : str
            任务类型
        max_show : int
            最多显示的失败任务数
        """
        failed_paths = results.get('failed', [])

        if not failed_paths:
            self.console.print(f"\n[green]没有失败的{task_type}任务！[/]")
            return

        self.console.print(f"\n[bold red]失败的{task_type}任务 (显示前{min(max_show, len(failed_paths))}个):[/]")

        table = Table(box=box.SIMPLE)
        table.add_column("序号", style="dim", width=6)
        table.add_column("路径", style="yellow")

        for i, path in enumerate(failed_paths[:max_show], 1):
            table.add_row(str(i), str(path.relative_to(self.base_dir)))

        self.console.print(table)

        if len(failed_paths) > max_show:
            self.console.print(f"[dim]... 还有 {len(failed_paths) - max_show} 个失败任务[/]")


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(description="检查QE任务完成状态")
    parser.add_argument("-d", "--depth", type=int, default=None,
                        help="最大搜索深度（默认：无限制）")
    parser.add_argument("-t", "--types", nargs='+', default=['relax', 'scf'],
                        help="要检查的任务类型（默认：relax scf）")
    parser.add_argument("--no-save", action='store_true',
                        help="不保存结果文件")
    parser.add_argument("--show-failed", type=int, default=10,
                        help="显示失败任务的数量（默认：10）")

    args = parser.parse_args()

    # 创建检查器
    checker = QETaskChecker(max_depth=args.depth)

    # 构建输出文件列表
    output_files = [f"{t}.out" for t in args.types]

    # 检查所有类型
    all_results = checker.check_all_types(output_files)

    # 打印和保存结果
    for task_type, results in all_results.items():
        checker.console.print("\n")
        checker.console.rule(f"[bold blue]{task_type.upper()}", style="blue")

        # 打印摘要
        checker.print_summary(results, task_type)

        # 打印失败详情
        if args.show_failed > 0:
            checker.print_failed_details(results, task_type, args.show_failed)

        # 保存结果
        if not args.no_save:
            checker.save_results(results, task_type)

    checker.console.print("\n[bold green]检查完成！[/]")


if __name__ == "__main__":
    main()
