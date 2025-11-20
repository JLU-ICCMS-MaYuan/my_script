"""
声子计算任务检查器

检查声子计算任务的完成情况，包括：
- 动力学矩阵文件（dyn）
- 电声耦合数据（elph_dir）
- 虚频检测
- 不可约表示进度

重构自: mytoolkit/qe_toolkit/cqepc.py

作者：Claude
创建时间：2025-11-20
"""

import os
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum

from rich.console import Console
from rich.table import Table
from rich import box


class DynStatus(str, Enum):
    """动力学矩阵状态"""
    DONE = "done"
    IMAGINARY = "imaginary"  # 有虚频
    INCOMPLETE = "incomplete"
    INEXISTENCE = "inexistence"


class ElphStatus(str, Enum):
    """电声耦合状态"""
    DONE = "done"
    EMPTY = "empty"
    INEXISTENCE = "inexistence"


@dataclass
class PhononTaskStatus:
    """单个q点的声子任务状态"""
    q_number: int
    dyn_status: DynStatus
    elph_status: ElphStatus
    representation_info: str = ""


class PhononTaskChecker:
    """声子计算任务检查器"""

    def __init__(self, work_dir: Path):
        """
        初始化检查器

        Parameters
        ----------
        work_dir : Path
            声子计算工作目录
        """
        self.work_dir = Path(work_dir)
        self.console = Console()
        self.prefix = None

    def _get_prefix(self) -> str:
        """获取dyn文件前缀"""
        if self.prefix:
            return self.prefix

        # 查找dyn0文件获取前缀
        dyn0_files = list(self.work_dir.glob("*.dyn0"))

        if not dyn0_files:
            raise FileNotFoundError("未找到dyn0文件，无法确定前缀")

        self.prefix = dyn0_files[0].stem.replace(".dyn0", "")
        return self.prefix

    def _check_dyn_file(self, q: int) -> Tuple[DynStatus, List[str]]:
        """
        检查动力学矩阵文件

        Parameters
        ----------
        q : int
            q点编号

        Returns
        -------
        status : DynStatus
            dyn文件状态
        freq_lines : List[str]
            频率行内容
        """
        q_dir = self.work_dir / str(q)

        if not q_dir.exists():
            return DynStatus.INEXISTENCE, []

        prefix = self._get_prefix()

        # 查找dyn文件（可能是dyn{q}或dyn）
        dyn_files = [
            q_dir / f"{prefix}.dyn{q}",
            q_dir / f"{prefix}.dyn",
        ]

        dyn_file = None
        for f in dyn_files:
            if f.exists() and f.stat().st_size > 0:
                dyn_file = f
                break

        if not dyn_file:
            return DynStatus.INEXISTENCE, []

        # 读取频率信息
        with open(dyn_file, 'r') as f:
            content = f.read()

        freq_lines = [line for line in content.split('\n') if 'freq' in line.lower()]

        if not freq_lines:
            return DynStatus.INCOMPLETE, []

        # 检查虚频
        has_imaginary = False
        for line in freq_lines:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    freq_value = float(parts[4])
                    if freq_value < 0:
                        has_imaginary = True
                        break
                except (ValueError, IndexError):
                    continue

        return (DynStatus.IMAGINARY if has_imaginary else DynStatus.DONE), freq_lines

    def _check_elph_dir(self, q: int) -> ElphStatus:
        """
        检查电声耦合目录

        Parameters
        ----------
        q : int
            q点编号

        Returns
        -------
        ElphStatus
            elph目录状态
        """
        q_dir = self.work_dir / str(q)
        elph_dir = q_dir / 'elph_dir'

        if not elph_dir.exists():
            return ElphStatus.INEXISTENCE

        if not any(elph_dir.iterdir()):
            return ElphStatus.EMPTY

        return ElphStatus.DONE

    def _get_representation_info(self, q: int) -> str:
        """
        获取不可约表示信息

        Parameters
        ----------
        q : int
            q点编号

        Returns
        -------
        str
            表示信息
        """
        q_dir = self.work_dir / str(q)
        ph_out = q_dir / 'split_ph.out'

        if not ph_out.exists():
            return ""

        # 查找最后一个Representation信息
        try:
            with open(ph_out, 'r') as f:
                lines = f.readlines()

            for line in reversed(lines):
                if "Representation #" in line:
                    return line.strip()
        except Exception:
            pass

        return ""

    def check_range(self, begin_q: int, end_q: int) -> List[PhononTaskStatus]:
        """
        检查q点范围内的任务状态

        Parameters
        ----------
        begin_q : int
            起始q点
        end_q : int
            结束q点

        Returns
        -------
        List[PhononTaskStatus]
            任务状态列表
        """
        results = []

        for q in range(begin_q, end_q + 1):
            # 检查dyn文件
            dyn_status, _ = self._check_dyn_file(q)

            # 检查elph目录
            elph_status = self._check_elph_dir(q)

            # 获取表示信息（如果dyn未完成）
            representation_info = ""
            if dyn_status in [DynStatus.INEXISTENCE, DynStatus.INCOMPLETE]:
                representation_info = self._get_representation_info(q)

            results.append(PhononTaskStatus(
                q_number=q,
                dyn_status=dyn_status,
                elph_status=elph_status,
                representation_info=representation_info
            ))

        return results

    def print_status(self, results: List[PhononTaskStatus]):
        """
        打印任务状态表格

        Parameters
        ----------
        results : List[PhononTaskStatus]
            任务状态列表
        """
        table = Table(
            title=f"声子任务状态 ({self.work_dir.name})",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold magenta",
        )

        table.add_column("Q点", style="cyan", width=10)
        table.add_column("DYN状态", width=20)
        table.add_column("ELPH状态", width=20)
        table.add_column("不可约表示", style="dim")

        for status in results:
            # DYN状态颜色
            if status.dyn_status == DynStatus.DONE:
                dyn_text = "[green]✓ 完成[/]"
            elif status.dyn_status == DynStatus.IMAGINARY:
                dyn_text = "[yellow]⚠ 虚频[/]"
            else:
                dyn_text = f"[red]✗ {status.dyn_status.value}[/]"

            # ELPH状态颜色
            if status.elph_status == ElphStatus.DONE:
                elph_text = "[green]✓ 完成[/]"
            elif status.elph_status == ElphStatus.EMPTY:
                elph_text = "[yellow]○ 空目录[/]"
            else:
                elph_text = "[dim]- 不存在[/]"

            table.add_row(
                str(status.q_number),
                dyn_text,
                elph_text,
                status.representation_info
            )

        self.console.print(table)

    def save_log(self, results: List[PhononTaskStatus], log_file: str = "phonon_check.log"):
        """
        保存检查日志

        Parameters
        ----------
        results : List[PhononTaskStatus]
            任务状态列表
        log_file : str
            日志文件名
        """
        log_path = self.work_dir / log_file

        with open(log_path, 'w') as f:
            f.write(f"{'Q点':<10} {'DYN状态':<20} {'ELPH状态':<20} {'不可约表示'}\n")
            f.write("-" * 80 + "\n")

            for status in results:
                f.write(
                    f"{status.q_number:<10} "
                    f"{status.dyn_status.value:<20} "
                    f"{status.elph_status.value:<20} "
                    f"{status.representation_info}\n"
                )

        self.console.print(f"\n日志已保存到: {log_path}")

    def get_summary(self, results: List[PhononTaskStatus]) -> Dict[str, int]:
        """
        获取统计摘要

        Parameters
        ----------
        results : List[PhononTaskStatus]
            任务状态列表

        Returns
        -------
        Dict[str, int]
            统计数据
        """
        summary = {
            'total': len(results),
            'dyn_done': 0,
            'dyn_imaginary': 0,
            'dyn_incomplete': 0,
            'elph_done': 0,
            'elph_empty': 0,
        }

        for status in results:
            if status.dyn_status == DynStatus.DONE:
                summary['dyn_done'] += 1
            elif status.dyn_status == DynStatus.IMAGINARY:
                summary['dyn_imaginary'] += 1
            else:
                summary['dyn_incomplete'] += 1

            if status.elph_status == ElphStatus.DONE:
                summary['elph_done'] += 1
            elif status.elph_status == ElphStatus.EMPTY:
                summary['elph_empty'] += 1

        return summary


def main():
    """命令行入口"""
    import sys

    if len(sys.argv) < 3:
        print("用法: python phonon_checker.py <begin_q> <end_q> [work_dir]")
        sys.exit(1)

    begin_q = int(sys.argv[1])
    end_q = int(sys.argv[2])
    work_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path.cwd()

    checker = PhononTaskChecker(work_dir)
    results = checker.check_range(begin_q, end_q)

    # 打印状态表格
    checker.print_status(results)

    # 保存日志
    checker.save_log(results)

    # 打印摘要
    summary = checker.get_summary(results)
    console = Console()
    console.print(f"\n[bold]摘要:[/]")
    console.print(f"  总q点数: {summary['total']}")
    console.print(f"  DYN完成: [green]{summary['dyn_done']}[/]")
    console.print(f"  DYN虚频: [yellow]{summary['dyn_imaginary']}[/]")
    console.print(f"  DYN未完成: [red]{summary['dyn_incomplete']}[/]")
    console.print(f"  ELPH完成: [green]{summary['elph_done']}[/]")


if __name__ == "__main__":
    main()
