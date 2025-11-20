"""
工作流基类模块

提供单步计算工作流的基础功能。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from abc import ABC, abstractmethod
from typing import Optional
import logging

from config.base import BaseConfig
from io.writers.base import BaseWriter

logger = logging.getLogger(__name__)


class BaseWorkflow(ABC):
    """
    工作流基类

    每个工作流代表一个完整的计算步骤（如relax、scf、phonon等）

    Attributes
    ----------
    config : BaseConfig
        配置对象
    work_dir : Path
        工作目录
    """

    def __init__(self, config: BaseConfig):
        """
        初始化工作流

        Parameters
        ----------
        config : BaseConfig
            配置对象
        """
        self.config = config
        self.work_dir = config.work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"初始化工作流: {self.__class__.__name__}")

    @abstractmethod
    def prepare_input(self) -> Path:
        """
        准备输入文件

        Returns
        -------
        Path
            输入文件路径
        """
        pass

    @abstractmethod
    def generate_submit_script(self) -> Path:
        """
        生成提交脚本

        Returns
        -------
        Path
            提交脚本路径
        """
        pass

    def run(self, submit: bool = True) -> bool:
        """
        执行工作流

        Parameters
        ----------
        submit : bool
            是否提交任务

        Returns
        -------
        bool
            是否成功
        """
        try:
            # 1. 准备输入文件
            input_file = self.prepare_input()
            logger.info(f"✓ 输入文件已生成: {input_file}")

            # 2. 生成提交脚本
            script_file = self.generate_submit_script()
            logger.info(f"✓ 提交脚本已生成: {script_file}")

            # 3. 提交任务（如果需要）
            if submit:
                self.submit_job(script_file)

            return True

        except Exception as e:
            logger.error(f"✗ 工作流执行失败: {e}")
            return False

    def submit_job(self, script_file: Path):
        """
        提交计算任务

        Parameters
        ----------
        script_file : Path
            提交脚本路径
        """
        # TODO: 实现任务提交逻辑
        logger.info(f"提交任务: {script_file}")
