"""
Pipeline基础架构

提供批量计算流程的核心功能。

作者：Claude
创建时间：2025-11-19
"""

from typing import List, Dict, Any, Optional
from pathlib import Path
from abc import ABC, abstractmethod
import logging

from core.types import TaskStatus
from scheduler.task import Task
from scheduler.executor import MixedParallelExecutor

logger = logging.getLogger(__name__)


class BasePipeline(ABC):
    """
    Pipeline基类

    用于定义批量计算流程。

    Attributes
    ----------
    structures : List[Path]
        结构文件列表
    work_dir : Path
        工作目录
    """

    def __init__(self, structures: List[Path], work_dir: Path, **config):
        """
        初始化Pipeline

        Parameters
        ----------
        structures : List[Path]
            结构文件路径列表
        work_dir : Path
            工作目录
        **config
            其他配置参数
        """
        self.structures = [Path(s) for s in structures]
        self.work_dir = Path(work_dir)
        self.config = config

        self.work_dir.mkdir(parents=True, exist_ok=True)

        # 任务列表
        self.tasks: List[Task] = []

        logger.info(f"初始化Pipeline: {len(structures)} 个结构")

    @abstractmethod
    def define_steps(self) -> List[Dict[str, Any]]:
        """
        定义流程步骤

        Returns
        -------
        List[Dict]
            步骤列表，每个步骤是一个字典：
            {
                'name': '步骤名称',
                'workflow': '工作流类名',
                'params': {...}
            }
        """
        pass

    def build_tasks(self):
        """构建任务列表"""
        steps = self.define_steps()

        for struct_file in self.structures:
            struct_dir = self.work_dir / struct_file.stem

            for step in steps:
                task = Task(
                    task_id=f"{struct_file.stem}_{step['name']}",
                    structure=struct_file,
                    step_name=step['name'],
                    work_dir=struct_dir,
                    config=step.get('params', {})
                )
                self.tasks.append(task)

        logger.info(f"已构建 {len(self.tasks)} 个任务")

    def run(self, max_workers: int = 4):
        """
        执行Pipeline

        Parameters
        ----------
        max_workers : int
            最大并行任务数
        """
        # 构建任务
        self.build_tasks()

        # 执行任务
        executor = MixedParallelExecutor(max_workers=max_workers)
        executor.execute(self.tasks)

        # 生成报告
        self.generate_report()

    def generate_report(self):
        """生成执行报告"""
        success_count = sum(1 for t in self.tasks if t.status == TaskStatus.SUCCESS)
        failed_count = sum(1 for t in self.tasks if t.status == TaskStatus.FAILED)

        report = f"""
Pipeline执行报告
{'='*60}
总任务数: {len(self.tasks)}
成功: {success_count}
失败: {failed_count}
成功率: {success_count/len(self.tasks)*100:.1f}%
{'='*60}
        """

        print(report)

        # 保存失败任务列表
        if failed_count > 0:
            failed_file = self.work_dir / 'failed_tasks.txt'
            with open(failed_file, 'w') as f:
                for task in self.tasks:
                    if task.status == TaskStatus.FAILED:
                        f.write(f"{task.task_id}: {task.error}\n")

            print(f"失败任务详情已保存到: {failed_file}")
