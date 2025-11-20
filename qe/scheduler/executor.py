"""
混合并行执行器

实现结构级并行 + 步骤内并行的任务执行。

作者：Claude
创建时间：2025-11-19
"""

from typing import List, Optional, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

from scheduler.task import Task
from core.types import TaskStatus

logger = logging.getLogger(__name__)


class MixedParallelExecutor:
    """
    混合并行执行器

    支持：
    1. 结构级并行：多个结构同时计算
    2. 步骤串行：每个结构的步骤按顺序执行
    """

    def __init__(self, max_workers: int = 4, progress_callback: Optional[Callable] = None):
        """
        初始化执行器

        Parameters
        ----------
        max_workers : int
            最大并行worker数量
        progress_callback : Callable, optional
            进度回调函数，参数为Task对象
        """
        self.max_workers = max_workers
        self.progress_callback = progress_callback

    def execute(self, tasks: List[Task]):
        """
        执行任务列表

        Parameters
        ----------
        tasks : List[Task]
            任务列表
        """
        logger.info(f"开始执行 {len(tasks)} 个任务，并行度={self.max_workers}")

        # 按结构分组任务
        task_groups = self._group_by_structure(tasks)

        # 并行执行每个结构的任务序列
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {}

            for struct_name, struct_tasks in task_groups.items():
                future = executor.submit(self._execute_structure_tasks, struct_tasks)
                futures[future] = struct_name

            # 等待所有任务完成
            for future in as_completed(futures):
                struct_name = futures[future]
                try:
                    future.result()
                    logger.info(f"✓ 结构 {struct_name} 完成")
                except Exception as e:
                    logger.error(f"✗ 结构 {struct_name} 失败: {e}")

    def _group_by_structure(self, tasks: List[Task]) -> dict:
        """按结构分组任务"""
        groups = {}
        for task in tasks:
            struct_name = task.structure.stem
            if struct_name not in groups:
                groups[struct_name] = []
            groups[struct_name].append(task)

        # 按步骤名称排序（确保步骤顺序）
        for tasks_list in groups.values():
            tasks_list.sort(key=lambda t: t.step_name)

        return groups

    def _execute_structure_tasks(self, tasks: List[Task]):
        """
        串行执行单个结构的所有步骤

        Parameters
        ----------
        tasks : List[Task]
            该结构的任务列表（已按步骤排序）
        """
        for task in tasks:
            try:
                self._execute_single_task(task)
            except Exception as e:
                # 失败则跳过该结构的后续步骤
                task.mark_failed(str(e))
                logger.error(f"任务失败，跳过后续步骤: {task.task_id}")
                break

    def _execute_single_task(self, task: Task):
        """
        执行单个任务

        Parameters
        ----------
        task : Task
            任务对象
        """
        task.mark_running()
        logger.info(f"⟳ 执行任务: {task.task_id}")

        # 通知进度监控
        if self.progress_callback:
            self.progress_callback(task)

        try:
            # TODO: 实际的任务执行逻辑
            # 这里需要调用相应的Workflow类
            # workflow = get_workflow(task.step_name)
            # workflow.run(task.work_dir, task.config)

            task.mark_success()
            logger.info(f"✓ 任务完成: {task.task_id}")

            # 通知进度监控
            if self.progress_callback:
                self.progress_callback(task)

        except Exception as e:
            task.mark_failed(str(e))
            logger.error(f"✗ 任务失败: {task.task_id} - {e}")

            # 通知进度监控
            if self.progress_callback:
                self.progress_callback(task)

            raise
