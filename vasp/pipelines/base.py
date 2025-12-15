"""
VASP Pipeline基类

提供Pipeline的基础功能：
- 任务状态管理
- 断点续传
- 日志记录
- 错误处理

作者：Claude
创建时间：2025-11-20
"""

import os
import json
import time
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Any
from enum import Enum

from vasp.utils.job import is_job_active

logger = logging.getLogger(__name__)


class StepStatus(str, Enum):
    """步骤状态枚举"""
    PENDING = "pending"       # 未执行
    PREPARED = "prepared"     # 仅生成输入
    RUNNING = "running"       # 正在执行
    COMPLETED = "completed"   # 已完成且有产物
    FAILED = "failed"         # 失败
    SKIPPED = "skipped"       # 跳过


class BasePipeline(ABC):
    """
    Pipeline基类

    所有VASP pipeline都应继承此类
    提供任务编排、状态管理、断点续传等核心功能
    """

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        checkpoint_file: Optional[Path] = None,
        report_file: Optional[Path] = None,
        max_retries: int = 1,
        retry_delay: int = 60,
        prepare_only: bool = False,
        pressure: float = 0.0,
    ):
        """
        初始化Pipeline

        Parameters
        ----------
        structure_file : Path
            输入结构文件（POSCAR格式）
        work_dir : Path
            工作目录
        checkpoint_file : Path, optional
            断点文件路径，默认为work_dir/pipeline_checkpoint.json
        max_retries : int
            任务失败时的最大重试次数
        retry_delay : int
            重试延迟（秒）
        """
        self.structure_file = Path(structure_file)
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        if self.structure_file.is_dir():
            raise ValueError(
                f"structure_file 必须是文件，当前是目录: {self.structure_file}。"
                " 可能是批量模式未正确识别，请检查输入目录下的结构后缀，或显式指定 --structure-ext。"
            )
        if not self.structure_file.exists():
            raise FileNotFoundError(f"结构文件不存在: {self.structure_file}")

        self.checkpoint_file = checkpoint_file or self.work_dir / "pipeline_checkpoint.json"
        self.report_file = report_file or self.work_dir / "pipeline_report.txt"
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.prepare_only = prepare_only
        self.pressure_gpa = pressure
        self.pressure_kbar = pressure * 10  # PSTRESS 单位 kBar

        # 步骤状态字典
        self.steps_status: Dict[str, StepStatus] = {}
        self.steps_data: Dict[str, Any] = {}

        # 加载断点（如果存在）
        self.log_file = self.work_dir / "pipeline.log"
        self._file_handler = None

        self._load_checkpoint()

    @abstractmethod
    def get_steps(self) -> List[str]:
        """
        返回Pipeline的所有步骤名称列表

        子类必须实现此方法
        """
        pass

    @abstractmethod
    def execute_step(self, step_name: str) -> bool:
        """
        执行单个步骤

        Parameters
        ----------
        step_name : str
            步骤名称

        Returns
        -------
        bool
            成功返回True，失败返回False
        """
        pass

    def run(self) -> bool:
        """
        执行整个Pipeline

        Returns
        -------
        bool
            所有步骤成功返回True，任意步骤失败返回False
        """
        # 每个任务单独日志文件
        root_logger = logging.getLogger()
        try:
            self._file_handler = logging.FileHandler(self.log_file, encoding="utf-8")
            self._file_handler.setLevel(logging.DEBUG)
            fmt = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            self._file_handler.setFormatter(fmt)
            root_logger.addHandler(self._file_handler)
        except Exception as e:  # pragma: no cover - 日志文件创建失败时仍继续
            logger.warning(f"创建日志文件失败: {e}")

        logger.info(f"开始执行Pipeline: {self.__class__.__name__}")
        logger.info(f"工作目录: {self.work_dir}")

        steps = self.get_steps()
        total_steps = len(steps)

        for idx, step_name in enumerate(steps, 1):
            logger.info(f"\n{'='*60}")
            logger.info(f"步骤 {idx}/{total_steps}: {step_name}")
            logger.info(f"{'='*60}")

            # 检查是否已完成
            if self.steps_status.get(step_name) == StepStatus.COMPLETED:
                logger.info(f"步骤 '{step_name}' 已完成，跳过")
                continue

            # 执行步骤（带重试机制）
            success = self._execute_step_with_retry(step_name)

            if not success:
                logger.error(f"Pipeline失败于步骤: {step_name}")
                self.steps_status[step_name] = StepStatus.FAILED
                self._save_checkpoint()
                return False

            # 标记完成并保存断点
            self.steps_status[step_name] = StepStatus.COMPLETED
            self._save_checkpoint()

            # 仅准备模式：生成输入但不提交，仅标记为 PREPARED
            if self.prepare_only:
                self.steps_status[step_name] = StepStatus.PREPARED
                self._save_checkpoint()
                logger.info("prepare_only=True，本次仅生成输入和脚本，不提交。继续准备后续步骤的输入。")
                continue

        logger.info(f"\n{'='*60}")
        logger.info("Pipeline执行完成！")
        logger.info(f"{'='*60}")

        # 生成最终报告
        self._generate_report()

        # 生成完成标记
        self._mark_finished()

        return True
    def __del__(self):
        # 清理文件日志句柄，避免重复记录
        if getattr(self, "_file_handler", None):
            root_logger = logging.getLogger()
            root_logger.removeHandler(self._file_handler)
            self._file_handler.close()

    def _execute_step_with_retry(self, step_name: str) -> bool:
        """
        执行步骤（带重试机制）

        Parameters
        ----------
        step_name : str
            步骤名称

        Returns
        -------
        bool
            成功返回True
        """
        for attempt in range(1, self.max_retries + 1):
            logger.info(f"尝试 {attempt}/{self.max_retries}")

            self.steps_status[step_name] = StepStatus.RUNNING

            try:
                success = self.execute_step(step_name)

                if success:
                    return True
                else:
                    logger.warning(f"步骤 '{step_name}' 执行失败")

            except Exception as e:
                logger.error(f"步骤 '{step_name}' 抛出异常: {e}", exc_info=True)

            # 如果不是最后一次尝试，等待后重试
            if attempt < self.max_retries:
                logger.info(f"等待 {self.retry_delay} 秒后重试...")
                time.sleep(self.retry_delay)

        return False

    def _wait_for_job(
        self,
        job_id: str,
        work_path: Path,
        queue_system: Optional[str] = None,
        check_interval: int = 30,
        timeout: int = 86400
    ) -> bool:
        """
        等待任务完成

        Parameters
        ----------
        job_id : str
            任务ID
        work_path : Path
            任务工作目录
        queue_system : str
            队列系统，便于查询作业状态
        check_interval : int
            检查间隔（秒）
        timeout : int
            超时时间（秒），默认24小时

        Returns
        -------
        bool
            任务成功完成返回True
        """
        logger.info(f"等待任务完成: {job_id}")

        start_time = time.time()

        while True:
            elapsed = time.time() - start_time

            # 检查超时
            if elapsed > timeout:
                logger.error(f"任务超时（{timeout}秒）")
                return False

            # 检查是否完成
            if self._check_job_completed(work_path):
                logger.info(f"任务完成！耗时: {elapsed:.0f}秒")
                return True

            # 队列状态检查（非bash）
            active = is_job_active(job_id, queue_system)
            if active is False:
                # 任务号消失，检查结果文件是否完成
                if self._check_job_completed(work_path):
                    logger.info(f"任务已结束且检测到完成标志，耗时: {elapsed:.0f}秒")
                    return True
                logger.error(f"队列中未找到任务 {job_id}，且未检测到完成标志，可能失败或被取消")
                return False

            # 等待
            time.sleep(check_interval)

    def _check_job_completed(self, work_path: Path) -> bool:
        """
        检查VASP任务是否完成

        Parameters
        ----------
        work_path : Path
            任务工作目录

        Returns
        -------
        bool
            完成返回True
        """
        # 检查OUTCAR文件中是否有"reached required accuracy"
        outcar = work_path / "OUTCAR"

        if not outcar.exists():
            return False

        try:
            with open(outcar, 'r') as f:
                content = f.read()

            # 检查是否正常结束
            if "reached required accuracy" in content or "writing wavefunctions" in content:
                return True

        except Exception as e:
            logger.warning(f"读取OUTCAR失败: {e}")

        return False

    def _check_convergence(self, work_path: Path) -> bool:
        """
        检查计算收敛性

        Parameters
        ----------
        work_path : Path
            任务工作目录

        Returns
        -------
        bool
            收敛返回True
        """
        outcar = work_path / "OUTCAR"

        if not outcar.exists():
            return False

        try:
            with open(outcar, 'r') as f:
                content = f.read()

            # 检查能量收敛
            if "reached required accuracy" in content:
                return True

            # 检查是否有WARNING
            if "WARNING" in content:
                logger.warning(f"OUTCAR中包含WARNING，请检查: {outcar}")

        except Exception as e:
            logger.error(f"读取OUTCAR失败: {e}")

        return False

    def _save_checkpoint(self):
        """保存断点状态"""
        checkpoint = {
            "steps_status": {k: v.value if isinstance(v, StepStatus) else v
                           for k, v in self.steps_status.items()},
            "steps_data": self.steps_data,
        }

        with open(self.checkpoint_file, 'w') as f:
            json.dump(checkpoint, f, indent=2)

        logger.debug(f"断点已保存: {self.checkpoint_file}")

    def _load_checkpoint(self):
        """加载断点状态"""
        if not self.checkpoint_file.exists():
            logger.info("未找到断点文件，从头开始")
            return

        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)

            self.steps_status = {k: StepStatus(v) for k, v in checkpoint.get("steps_status", {}).items()}
            self.steps_data = checkpoint.get("steps_data", {})

            logger.info(f"已加载断点: {self.checkpoint_file}")
            logger.info(f"已完成步骤: {[k for k, v in self.steps_status.items() if v == StepStatus.COMPLETED]}")

        except Exception as e:
            logger.warning(f"加载断点失败: {e}，从头开始")

    def _generate_report(self):
        """生成Pipeline执行报告"""
        with open(self.report_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write(f"Pipeline执行报告: {self.__class__.__name__}\n")
            f.write("="*60 + "\n\n")

            f.write(f"工作目录: {self.work_dir}\n")
            f.write(f"输入结构: {self.structure_file}\n\n")

            f.write("步骤执行状态:\n")
            f.write("-"*60 + "\n")
            for step, status in self.steps_status.items():
                f.write(f"  {step:30s} : {status.value}\n")

            f.write("\n" + "="*60 + "\n")

        logger.info(f"报告已生成: {self.report_file}")

    def _mark_finished(self):
        """在工作目录写入完成标记文件。"""
        try:
            marker = self.work_dir / "finished"
            marker.write_text("finished\n")
        except Exception as exc:  # pragma: no cover
            logger.warning(f"写入完成标记失败: {exc}")
