"""
统一日志模块

适用于QE、VASP、EPW三个计算软件的通用日志配置。
消除原有三个项目90%的重复代码。

作者：Claude
创建时间：2025-11-20
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Union


# 日志级别映射
LOG_LEVELS = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL,
}


def setup_logger(
    name: str,
    level: Union[str, int] = 'INFO',
    log_file: Optional[Path] = None,
    console: bool = True,
    format_string: Optional[str] = None
) -> logging.Logger:
    """
    配置统一的日志记录器

    Parameters
    ----------
    name : str
        Logger名称
    level : str or int
        日志级别 ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')
    log_file : Path, optional
        日志文件路径
    console : bool
        是否输出到控制台
    format_string : str, optional
        自定义日志格式

    Returns
    -------
    logging.Logger
        配置好的logger对象

    Examples
    --------
    >>> logger = setup_logger('qe', 'DEBUG', log_file=Path('qe.log'))
    >>> logger.info('开始计算')
    """
    # 获取logger
    logger = logging.getLogger(name)
    logger.setLevel(LOG_LEVELS.get(level, logging.INFO) if isinstance(level, str) else level)

    # 清除已有的handlers（避免重复）
    logger.handlers.clear()

    # 日志格式
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    formatter = logging.Formatter(format_string)

    # 控制台handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logger.level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # 文件handler
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logger.level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # 防止日志传播到root logger
    logger.propagate = False

    return logger


def setup_default_logger(software: str, level: str = 'INFO', work_dir: Optional[Path] = None) -> logging.Logger:
    """
    为特定软件设置默认logger

    Parameters
    ----------
    software : str
        软件名称 ('qe', 'vasp', 'epw')
    level : str
        日志级别
    work_dir : Path, optional
        工作目录（日志文件将保存在此目录）

    Returns
    -------
    logging.Logger
        配置好的logger

    Examples
    --------
    >>> logger = setup_default_logger('qe', 'DEBUG', Path('./calculations'))
    """
    log_file = None
    if work_dir:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        log_file = work_dir / f"{software}.log"

    return setup_logger(
        name=software,
        level=level,
        log_file=log_file,
        console=True
    )


class LoggerAdapter:
    """
    Logger适配器，用于向后兼容旧代码

    将原有的logging_level参数转换为新的logger。

    Examples
    --------
    >>> adapter = LoggerAdapter('qe')
    >>> logger = adapter.get_logger('INFO', Path('./work'))
    """

    def __init__(self, software: str):
        """
        初始化适配器

        Parameters
        ----------
        software : str
            软件名称
        """
        self.software = software

    def get_logger(self, logging_level: str = 'INFO', work_path: Optional[Path] = None) -> logging.Logger:
        """
        获取logger（兼容旧接口）

        Parameters
        ----------
        logging_level : str
            日志级别
        work_path : Path, optional
            工作目录

        Returns
        -------
        logging.Logger
            Logger对象
        """
        return setup_default_logger(self.software, logging_level, work_path)


# 便捷函数
def get_qe_logger(level: str = 'INFO', work_dir: Optional[Path] = None) -> logging.Logger:
    """获取QE logger"""
    return setup_default_logger('qe', level, work_dir)


def get_vasp_logger(level: str = 'INFO', work_dir: Optional[Path] = None) -> logging.Logger:
    """获取VASP logger"""
    return setup_default_logger('vasp', level, work_dir)


def get_epw_logger(level: str = 'INFO', work_dir: Optional[Path] = None) -> logging.Logger:
    """获取EPW logger"""
    return setup_default_logger('epw', level, work_dir)
