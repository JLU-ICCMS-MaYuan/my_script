"""
配置基类模块

所有QE配置类的基类，提供参数验证、默认值管理、配置导出等通用功能。

作者：Claude
创建时间：2025-11-19
"""

from typing import Dict, Any, Optional, List
from pathlib import Path
from abc import ABC, abstractmethod
import logging

from core.types import PathLike, to_path
from core.exceptions import InvalidParameterError, MissingParameterError

logger = logging.getLogger(__name__)


class BaseConfig(ABC):
    """
    QE配置基类

    所有配置类都应该继承此类，提供统一的配置管理接口。
    """

    def __init__(self, work_dir: PathLike, input_file: PathLike, pressure: float = 0.0):
        """初始化基础配置"""
        self.work_dir = to_path(work_dir)
        self.input_file = to_path(input_file)
        self.pressure = pressure

        self.work_dir.mkdir(parents=True, exist_ok=True)
        self._params: Dict[str, Any] = {}
        self._set_defaults()

        logger.info(f"初始化配置: {self.__class__.__name__}")

    @abstractmethod
    def _set_defaults(self):
        """设置默认参数值（子类实现）"""
        pass

    @abstractmethod
    def validate(self) -> bool:
        """验证配置参数（子类实现）"""
        pass

    def set_param(self, param: str, value: Any):
        """设置参数值"""
        setattr(self, param, value)
        self._params[param] = value

    def get_param(self, param: str, default: Any = None) -> Any:
        """获取参数值"""
        return getattr(self, param, default)

    def to_dict(self) -> Dict[str, Any]:
        """导出为字典"""
        return {
            'work_dir': str(self.work_dir),
            'input_file': str(self.input_file),
            'pressure': self.pressure,
            **self._params
        }

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "BaseConfig":
        """从字典创建配置"""
        work_dir = config_dict.pop('work_dir')
        input_file = config_dict.pop('input_file')
        pressure = config_dict.pop('pressure', 0.0)

        instance = cls(work_dir, input_file, pressure)
        for key, value in config_dict.items():
            instance.set_param(key, value)

        return instance


# 默认参数配置
DEFAULT_PARAMS = {
    'ecutwfc': 80.0,
    'ecutrho': 320.0,
    'conv_thr': '1.0d-8',
    'mixing_beta': 0.3,
    'degauss': 0.05,
}
