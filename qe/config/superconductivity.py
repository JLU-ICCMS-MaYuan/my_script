"""
超导性质计算配置模块

作者：Claude
创建时间：2025-11-19
"""

from typing import List
from config.base import BaseConfig


class SuperconductivityConfig(BaseConfig):
    """超导性质计算配置"""

    def _set_defaults(self):
        """设置默认参数"""
        # 屏蔽常数列表
        self.set_param('mu_values', [0.10, 0.13])

        # Eliashberg求解参数
        self.set_param('temperature_steps', 100)  # 温度点数

        # Lambda计算参数
        self.set_param('top_freq', 0.0)  # 最高频率（0表示自动）
        self.set_param('broaden', 0.5)   # 展宽方法
        self.set_param('smearing_method', 1)

        # 数据来源
        self.set_param('a2fdos', True)   # 从a2F.dos读取
        self.set_param('alpha2fdat', False)  # 从alpha2F.dat读取

    def validate(self) -> bool:
        """验证参数"""
        mu_values = self.get_param('mu_values')
        if not mu_values or not isinstance(mu_values, list):
            from core.exceptions import InvalidParameterError
            raise InvalidParameterError('mu_values', mu_values, '必须是列表')

        return True
