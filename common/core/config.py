"""
统一配置读取模块

适用于QE、VASP、EPW三个计算软件的通用配置读取。
消除原有三个项目95%的重复代码。

作者：Claude
创建时间：2025-11-20
"""

from argparse import ArgumentParser
from typing import Dict, Any, Optional
from abc import ABC, abstractmethod


class BaseConfig(ABC):
    """
    配置读取基类

    统一处理命令行参数到配置字典的转换。
    子类只需指定workflow_key即可。

    Examples
    --------
    >>> class QEConfig(BaseConfig):
    ...     workflow_key = "qe_workflow"
    >>>
    >>> args = parser.parse_args()
    >>> config = QEConfig(args).read_config()
    """

    # 子类需要指定workflow参数的key名
    workflow_key: str = None

    def __init__(self, args: ArgumentParser):
        """
        初始化配置读取器

        Parameters
        ----------
        args : ArgumentParser
            命令行参数对象
        """
        self.args = args

    def read_config(self) -> Dict[str, Any]:
        """
        读取配置

        Returns
        -------
        Dict[str, Any]
            配置字典
        """
        config = {}

        # 通用参数
        config["input_file_path"] = self.args.input_file_path
        config["work_path"] = self.args.work_path
        config["submit_job_system"] = self.args.submit_job_system
        config["logging_level"] = self.args.logging_level

        # 项目特定参数
        self._add_project_specific(config)

        # Workflow参数（子类指定）
        if self.workflow_key:
            config["workflow"] = getattr(self.args, self.workflow_key, None)

        # 额外参数（key=value格式）
        if hasattr(self.args, 'more_args') and self.args.more_args:
            for other_arg in self.args.more_args:
                arg_name, value = other_arg.split("=", 1)
                config[arg_name] = value

        return config

    def _add_project_specific(self, config: Dict[str, Any]):
        """
        添加项目特定参数

        子类可重写此方法添加特定参数。

        Parameters
        ----------
        config : Dict[str, Any]
            配置字典
        """
        pass


class QEConfig(BaseConfig):
    """QE配置读取器"""
    workflow_key = "qe_workflow"

    def _add_project_specific(self, config: Dict[str, Any]):
        """添加QE特定参数"""
        if hasattr(self.args, 'press'):
            config["press"] = self.args.press
        if hasattr(self.args, 'pp_dir'):
            config["pp_dir"] = self.args.pp_dir


class VASPConfig(BaseConfig):
    """VASP配置读取器"""
    workflow_key = "vasp_workflow"

    def _add_project_specific(self, config: Dict[str, Any]):
        """添加VASP特定参数"""
        if hasattr(self.args, 'press'):
            config["press"] = self.args.press
        if hasattr(self.args, 'presses'):
            config["presses"] = self.args.presses
        if hasattr(self.args, 'pp_dir'):
            config["pp_dir"] = self.args.pp_dir


class EPWConfig(BaseConfig):
    """EPW配置读取器"""
    workflow_key = "epw_workflow"


# 便捷函数
def create_config(software: str, args: ArgumentParser) -> Dict[str, Any]:
    """
    根据软件类型创建配置

    Parameters
    ----------
    software : str
        软件类型: 'qe', 'vasp', 'epw'
    args : ArgumentParser
        命令行参数

    Returns
    -------
    Dict[str, Any]
        配置字典

    Examples
    --------
    >>> config = create_config('qe', args)
    """
    config_classes = {
        'qe': QEConfig,
        'vasp': VASPConfig,
        'epw': EPWConfig,
    }

    config_class = config_classes.get(software.lower())
    if not config_class:
        raise ValueError(f"不支持的软件类型: {software}")

    return config_class(args).read_config()
