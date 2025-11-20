"""
统一任务提交模块

适用于QE、VASP、EPW三个计算软件的通用任务提交。
消除原有三个项目85%的重复代码。

支持的队列系统：
- SLURM (sbatch)
- PBS (qsub)
- LSF (bsub)
- Bash (后台执行)

作者：Claude
创建时间：2025-11-20
"""

import os
import re
import logging
from pathlib import Path
from typing import List, Optional, Union
from enum import Enum


logger = logging.getLogger(__name__)


class QueueSystem(str, Enum):
    """队列系统类型"""
    SLURM = "slurm"
    PBS = "pbs"
    LSF = "lsf"
    BASH = "bash"


class JobSubmitter:
    """
    统一的任务提交器

    封装不同队列系统的任务提交逻辑。

    Parameters
    ----------
    queue_system : str or QueueSystem
        队列系统类型
    work_path : Path, optional
        工作目录

    Examples
    --------
    >>> submitter = JobSubmitter("slurm", work_path=Path("./calc"))
    >>> jobids = submitter.submit_job("submit.sh")
    >>> print(f"Job IDs: {jobids}")
    """

    # 队列系统对应的提交命令
    SUBMIT_COMMANDS = {
        QueueSystem.SLURM: "sbatch",
        QueueSystem.PBS: "qsub",
        QueueSystem.LSF: "bsub <",
        QueueSystem.BASH: "bash",
    }

    def __init__(self, queue_system: Union[str, QueueSystem], work_path: Optional[Path] = None):
        """初始化任务提交器"""
        if isinstance(queue_system, str):
            queue_system = QueueSystem(queue_system.lower())

        self.queue_system = queue_system
        self.work_path = Path(work_path) if work_path else Path.cwd()

        # 获取提交命令
        self.submit_order = self.SUBMIT_COMMANDS.get(queue_system, '')

        if not self.submit_order:
            raise ValueError(f"不支持的队列系统: {queue_system}")

    def submit_job(
        self,
        job_script: Union[str, Path],
        work_dir: Optional[Path] = None,
        required_files: Optional[List[str]] = None
    ) -> List[str]:
        """
        提交任务

        Parameters
        ----------
        job_script : str or Path
            任务脚本文件名
        work_dir : Path, optional
            工作目录（默认使用初始化时的work_path）
        required_files : List[str], optional
            必需的文件列表，提交前检查

        Returns
        -------
        List[str]
            Job IDs（任务ID或进程ID列表）

        Raises
        ------
        FileNotFoundError
            如果必需文件不存在
        RuntimeError
            如果任务提交失败
        """
        work_dir = Path(work_dir) if work_dir else self.work_path

        # 检查任务脚本是否存在
        job_file = work_dir / job_script
        if not job_file.exists():
            raise FileNotFoundError(f"任务脚本不存在: {job_file}")

        # 检查必需文件
        if required_files:
            for required_file in required_files:
                file_path = work_dir / required_file
                if not file_path.exists():
                    raise FileNotFoundError(f"必需文件不存在: {file_path}")

        # 切换到工作目录
        cwd = os.getcwd()
        os.chdir(work_dir)

        try:
            # 提交任务
            if self.queue_system == QueueSystem.BASH:
                # Bash后台执行
                logger.debug(f"nohup {self.submit_order} {job_script} > bash.log 2>&1 &")
                os.popen(f"nohup {self.submit_order} {job_script} > bash.log 2>&1 &").read()
                jobids = self._get_pid()
            else:
                # 队列系统提交
                logger.debug(f"{self.submit_order} {job_script}")
                res = os.popen(f"{self.submit_order} {job_script}").read()
                jobids = re.findall(r"\d+", res)

            logger.info(f"{job_script} 已提交。Job IDs: {jobids}")

            # 验证提交成功
            if not jobids:
                raise RuntimeError(
                    f"任务提交失败！Job IDs为空。\n"
                    f"队列系统: {self.queue_system}\n"
                    f"提交命令: {self.submit_order} {job_script}"
                )

            return jobids

        finally:
            # 恢复原目录
            os.chdir(cwd)

    def submit_direct_command(
        self,
        command: str,
        work_dir: Optional[Path] = None,
        background: bool = False
    ) -> Optional[List[str]]:
        """
        直接执行命令（不通过队列）

        Parameters
        ----------
        command : str
            要执行的命令
        work_dir : Path, optional
            工作目录
        background : bool
            是否后台执行

        Returns
        -------
        List[str] or None
            如果后台执行返回PIDs，否则返回None
        """
        work_dir = Path(work_dir) if work_dir else self.work_path

        cwd = os.getcwd()
        os.chdir(work_dir)

        try:
            if background:
                logger.debug(f"nohup {command} > command.log 2>&1 &")
                os.popen(f"nohup {command} > command.log 2>&1 &").read()
                return self._get_pid()
            else:
                logger.debug(f"{command}")
                os.system(command)
                return None

        finally:
            os.chdir(cwd)

    def _get_pid(self) -> List[str]:
        """
        获取当前bash后台任务的PID

        Returns
        -------
        List[str]
            PID列表
        """
        import subprocess

        try:
            # 使用ps命令获取当前用户的bash进程
            result = subprocess.run(
                ['ps', '-u', os.getenv('USER'), '-o', 'pid,cmd'],
                capture_output=True,
                text=True
            )

            # 从输出中提取bash相关的PID
            lines = result.stdout.strip().split('\n')[1:]  # 跳过表头
            pids = []

            for line in lines:
                parts = line.strip().split(None, 1)
                if len(parts) == 2:
                    pid, cmd = parts
                    if 'bash' in cmd.lower():
                        pids.append(pid)

            return pids if pids else ['unknown']

        except Exception as e:
            logger.warning(f"获取PID失败: {e}")
            return ['unknown']

    def kill_job(self, jobids: Union[str, List[str]]):
        """
        终止任务

        Parameters
        ----------
        jobids : str or List[str]
            要终止的任务ID
        """
        if isinstance(jobids, str):
            jobids = [jobids]

        for jobid in jobids:
            if self.queue_system == QueueSystem.SLURM:
                os.system(f"scancel {jobid}")
            elif self.queue_system == QueueSystem.PBS:
                os.system(f"qdel {jobid}")
            elif self.queue_system == QueueSystem.LSF:
                os.system(f"bkill {jobid}")
            elif self.queue_system == QueueSystem.BASH:
                os.system(f"kill -9 {jobid}")

            logger.info(f"已终止任务: {jobid}")

    def check_job_status(self, jobid: str) -> bool:
        """
        检查任务是否仍在运行

        Parameters
        ----------
        jobid : str
            任务ID

        Returns
        -------
        bool
            是否仍在运行
        """
        if self.queue_system == QueueSystem.SLURM:
            result = os.popen(f"squeue -j {jobid}").read()
            return jobid in result
        elif self.queue_system == QueueSystem.PBS:
            result = os.popen(f"qstat {jobid}").read()
            return jobid in result
        elif self.queue_system == QueueSystem.LSF:
            result = os.popen(f"bjobs {jobid}").read()
            return jobid in result
        elif self.queue_system == QueueSystem.BASH:
            result = os.popen(f"ps -p {jobid}").read()
            return jobid in result

        return False


# 便捷函数
def submit_job_simple(
    queue_system: str,
    job_script: str,
    work_dir: Path,
    **kwargs
) -> List[str]:
    """
    简化的任务提交函数

    Parameters
    ----------
    queue_system : str
        队列系统类型
    job_script : str
        任务脚本
    work_dir : Path
        工作目录
    **kwargs
        传递给JobSubmitter.submit_job的其他参数

    Returns
    -------
    List[str]
        Job IDs

    Examples
    --------
    >>> jobids = submit_job_simple("slurm", "submit.sh", Path("./calc"))
    """
    submitter = JobSubmitter(queue_system, work_dir)
    return submitter.submit_job(job_script, **kwargs)
