"""
任务定义模块

作者：Claude
创建时间：2025-11-19
"""

from typing import Optional, List, Dict, Any
from pathlib import Path
from dataclasses import dataclass, field
from datetime import datetime

from core.types import TaskStatus


@dataclass
class Task:
    """
    单个计算任务

    Attributes
    ----------
    task_id : str
        任务唯一标识
    structure : Path
        结构文件路径
    step_name : str
        步骤名称
    work_dir : Path
        工作目录
    config : Dict
        配置参数
    """
    task_id: str
    structure: Path
    step_name: str
    work_dir: Path
    config: Dict[str, Any] = field(default_factory=dict)

    # 任务状态
    status: TaskStatus = TaskStatus.PENDING
    dependencies: List[str] = field(default_factory=list)

    # 执行信息
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    error: Optional[str] = None

    def mark_running(self):
        """标记为运行中"""
        self.status = TaskStatus.RUNNING
        self.start_time = datetime.now()

    def mark_success(self):
        """标记为成功"""
        self.status = TaskStatus.SUCCESS
        self.end_time = datetime.now()

    def mark_failed(self, error: str):
        """标记为失败"""
        self.status = TaskStatus.FAILED
        self.end_time = datetime.now()
        self.error = error

    @property
    def duration(self) -> Optional[float]:
        """执行时长（秒）"""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return None

    def __repr__(self) -> str:
        return f"Task({self.task_id}, {self.status.value})"
