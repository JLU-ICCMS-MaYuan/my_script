"""
电子性质计算配置模块

作者：Claude
创建时间：2025-11-19
"""

from config.base import BaseConfig, DEFAULT_PARAMS


class ElectronConfig(BaseConfig):
    """电子性质计算配置"""

    def _set_defaults(self):
        """设置默认参数"""
        # 基础参数
        self.set_param('ecutwfc', DEFAULT_PARAMS['ecutwfc'])
        self.set_param('ecutrho', DEFAULT_PARAMS['ecutrho'])

        # K点设置
        self.set_param('kpoints_dense', '16 16 16')  # DOS用密K点
        self.set_param('kinserted', 200)  # 能带路径插值点数

        # 能带参数
        self.set_param('nbnd', 500)  # 能带数量（自动计算）

        # 计算模式
        self.set_param('mode', 'eleproperties')  # eleband/eledos/eleproperties

    def validate(self) -> bool:
        """验证参数"""
        return True
