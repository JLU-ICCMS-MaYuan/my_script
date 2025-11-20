"""
实时进度监控模块

使用rich库提供美观的命令行进度显示，包括：
- 进度条
- 任务状态表格
- 实时日志
- 统计信息

作者：Claude
创建时间：2025-11-20
"""

from typing import Dict, List, Optional
from pathlib import Path
from datetime import datetime, timedelta
import time

from rich.console import Console
from rich.table import Table
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeRemainingColumn,
    TimeElapsedColumn,
)
from rich.live import Live
from rich.layout import Layout
from rich.panel import Panel
from rich.text import Text
from rich import box

from scheduler.task import Task, TaskStatus


class ProgressMonitor:
    """
    实时进度监控器

    功能：
    - 显示整体进度条
    - 显示任务状态表格
    - 实时更新计算进度
    - 统计成功/失败/运行中任务
    """

    def __init__(self, total_tasks: int, refresh_rate: float = 0.5):
        """
        初始化进度监控器

        Parameters
        ----------
        total_tasks : int
            总任务数
        refresh_rate : float
            刷新率（秒）
        """
        self.console = Console()
        self.total_tasks = total_tasks
        self.refresh_rate = refresh_rate
        self.start_time = datetime.now()

        # 任务统计
        self.stats = {
            'pending': 0,
            'running': 0,
            'success': 0,
            'failed': 0
        }

        # 任务列表
        self.tasks: List[Task] = []

        # Rich组件
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(bar_width=40),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=self.console,
        )

        self.main_task_id = None

    def add_tasks(self, tasks: List[Task]):
        """添加要监控的任务"""
        self.tasks = tasks
        self._update_stats()

    def _update_stats(self):
        """更新统计信息"""
        self.stats = {
            'pending': 0,
            'running': 0,
            'success': 0,
            'failed': 0
        }

        for task in self.tasks:
            status = task.status.value
            if status in self.stats:
                self.stats[status] += 1

    def _create_status_table(self) -> Table:
        """创建任务状态表格"""
        table = Table(
            title="任务状态",
            box=box.ROUNDED,
            show_header=True,
            header_style="bold magenta",
            expand=False,
        )

        table.add_column("结构", style="cyan", width=20)
        table.add_column("步骤", style="yellow", width=15)
        table.add_column("状态", width=12)
        table.add_column("进度", width=10)
        table.add_column("耗时", width=10)

        # 只显示最近的20个任务
        recent_tasks = self.tasks[-20:] if len(self.tasks) > 20 else self.tasks

        for task in recent_tasks:
            # 状态颜色
            if task.status == TaskStatus.SUCCESS:
                status_text = Text("✓ 成功", style="green")
            elif task.status == TaskStatus.FAILED:
                status_text = Text("✗ 失败", style="red")
            elif task.status == TaskStatus.RUNNING:
                status_text = Text("▶ 运行中", style="yellow")
            else:
                status_text = Text("○ 等待", style="dim")

            # 进度百分比
            if task.status == TaskStatus.SUCCESS:
                progress = "100%"
            elif task.status == TaskStatus.RUNNING:
                progress = f"{task.progress:.0f}%"
            else:
                progress = "0%"

            # 计算耗时
            if task.start_time:
                if task.end_time:
                    elapsed = task.end_time - task.start_time
                else:
                    elapsed = datetime.now() - task.start_time
                elapsed_str = str(timedelta(seconds=int(elapsed.total_seconds())))
            else:
                elapsed_str = "-"

            table.add_row(
                task.structure_name[:18],
                task.step_name,
                status_text,
                progress,
                elapsed_str,
            )

        return table

    def _create_stats_panel(self) -> Panel:
        """创建统计信息面板"""
        completed = self.stats['success'] + self.stats['failed']
        total_progress = (completed / self.total_tasks * 100) if self.total_tasks > 0 else 0

        elapsed = datetime.now() - self.start_time
        elapsed_str = str(timedelta(seconds=int(elapsed.total_seconds())))

        # 估算剩余时间
        if completed > 0:
            avg_time_per_task = elapsed.total_seconds() / completed
            remaining_tasks = self.total_tasks - completed
            eta_seconds = avg_time_per_task * remaining_tasks
            eta_str = str(timedelta(seconds=int(eta_seconds)))
        else:
            eta_str = "计算中..."

        stats_text = f"""
[bold cyan]总任务数:[/] {self.total_tasks}
[bold green]已完成:[/] {self.stats['success']}
[bold red]失败:[/] {self.stats['failed']}
[bold yellow]运行中:[/] {self.stats['running']}
[bold dim]等待中:[/] {self.stats['pending']}

[bold magenta]总体进度:[/] {total_progress:.1f}%
[bold blue]已用时间:[/] {elapsed_str}
[bold blue]预计剩余:[/] {eta_str}
"""
        return Panel(stats_text, title="统计信息", border_style="cyan", box=box.ROUNDED)

    def _create_layout(self) -> Layout:
        """创建整体布局"""
        layout = Layout()

        layout.split_column(
            Layout(name="header", size=3),
            Layout(name="main"),
            Layout(name="footer", size=10),
        )

        layout["main"].split_row(
            Layout(name="tasks", ratio=2),
            Layout(name="stats", ratio=1),
        )

        # 设置各部分内容
        layout["header"].update(
            Panel(
                Text("QE计算任务监控", justify="center", style="bold white"),
                style="bold blue",
            )
        )

        layout["tasks"].update(self._create_status_table())
        layout["stats"].update(self._create_stats_panel())

        # 底部进度条
        progress_text = f"总体进度: {self.stats['success']}/{self.total_tasks} 完成"
        layout["footer"].update(
            Panel(
                Text(progress_text, justify="center"),
                border_style="green",
            )
        )

        return layout

    def update_task_status(self, task: Task):
        """
        更新单个任务状态

        Parameters
        ----------
        task : Task
            更新的任务
        """
        # 找到并更新任务
        for i, t in enumerate(self.tasks):
            if t.task_id == task.task_id:
                self.tasks[i] = task
                break

        self._update_stats()

    def run_with_monitor(self, executor_func, *args, **kwargs):
        """
        运行任务并显示监控界面

        Parameters
        ----------
        executor_func : callable
            执行器函数
        *args, **kwargs
            传递给执行器的参数
        """
        with Live(
            self._create_layout(),
            console=self.console,
            refresh_per_second=1 / self.refresh_rate,
            screen=True,
        ) as live:
            # 定义更新回调
            def update_callback(task: Task):
                self.update_task_status(task)
                live.update(self._create_layout())

            # 执行任务（传入回调）
            if 'progress_callback' in kwargs:
                original_callback = kwargs['progress_callback']

                def combined_callback(task):
                    update_callback(task)
                    original_callback(task)

                kwargs['progress_callback'] = combined_callback
            else:
                kwargs['progress_callback'] = update_callback

            # 运行执行器
            result = executor_func(*args, **kwargs)

            # 最后更新一次
            live.update(self._create_layout())
            time.sleep(1)  # 让用户看到最终结果

            return result

    def print_summary(self):
        """打印最终摘要"""
        self.console.print("\n")
        self.console.rule("[bold blue]计算完成", style="blue")

        # 成功率
        success_rate = (self.stats['success'] / self.total_tasks * 100) if self.total_tasks > 0 else 0

        summary_table = Table(
            title="最终统计",
            box=box.DOUBLE,
            show_header=True,
            header_style="bold cyan",
        )

        summary_table.add_column("项目", style="cyan")
        summary_table.add_column("数值", style="yellow", justify="right")

        summary_table.add_row("总任务数", str(self.total_tasks))
        summary_table.add_row("成功", f"[green]{self.stats['success']}[/]")
        summary_table.add_row("失败", f"[red]{self.stats['failed']}[/]")
        summary_table.add_row("成功率", f"[bold green]{success_rate:.1f}%[/]")

        elapsed = datetime.now() - self.start_time
        elapsed_str = str(timedelta(seconds=int(elapsed.total_seconds())))
        summary_table.add_row("总耗时", elapsed_str)

        self.console.print(summary_table)

        # 如果有失败任务，列出来
        if self.stats['failed'] > 0:
            self.console.print("\n[bold red]失败任务列表:[/]")

            failed_table = Table(box=box.SIMPLE)
            failed_table.add_column("结构", style="cyan")
            failed_table.add_column("步骤", style="yellow")
            failed_table.add_column("错误信息", style="red")

            for task in self.tasks:
                if task.status == TaskStatus.FAILED:
                    error_msg = task.error_message if task.error_message else "未知错误"
                    failed_table.add_row(
                        task.structure_name,
                        task.step_name,
                        error_msg[:50],  # 截断错误信息
                    )

            self.console.print(failed_table)


class SimpleProgressBar:
    """
    简单进度条（用于单个任务）

    适用于单个长时间运行的任务
    """

    def __init__(self, total: int, description: str = "处理中"):
        """
        初始化简单进度条

        Parameters
        ----------
        total : int
            总步数
        description : str
            描述文字
        """
        self.console = Console()
        self.total = total
        self.description = description

        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.description}"),
            BarColumn(bar_width=50),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            console=self.console,
        )

        self.task_id = None

    def __enter__(self):
        """进入上下文"""
        self.progress.start()
        self.task_id = self.progress.add_task(self.description, total=self.total)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """退出上下文"""
        self.progress.stop()

    def update(self, advance: int = 1, description: Optional[str] = None):
        """
        更新进度

        Parameters
        ----------
        advance : int
            前进步数
        description : str, optional
            更新描述
        """
        if description:
            self.progress.update(self.task_id, advance=advance, description=description)
        else:
            self.progress.update(self.task_id, advance=advance)

    def set_total(self, total: int):
        """设置新的总步数"""
        self.progress.update(self.task_id, total=total)
