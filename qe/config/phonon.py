"""
声子计算配置模块

作者：Claude
创建时间：2025-11-19
"""

from config.base import BaseConfig, DEFAULT_PARAMS


class PhononConfig(BaseConfig):
    """声子计算配置"""

    def _set_defaults(self):
        """设置默认参数"""
        # 基础参数
        self.set_param('ecutwfc', DEFAULT_PARAMS['ecutwfc'])
        self.set_param('ecutrho', DEFAULT_PARAMS['ecutrho'])

        # Q点设置
        self.set_param('qpoints', '6 6 6')

        # 声子计算参数
        self.set_param('tr2_ph', '1.0d-14')  # 声子收敛阈值
        self.set_param('alpha_mix', 0.7)     # 混合参数

        # 电声耦合参数
        self.set_param('el_ph_sigma', 0.005)  # 展宽起始值
        self.set_param('el_ph_nsigma', 10)    # 展宽数量

        # 计算模式
        self.set_param('mode', 'nosplit')     # nosplit/split_dyn0/split_assignQ

    def validate(self) -> bool:
        """验证参数"""
        # 验证Q点设置
        qpoints = self.get_param('qpoints')
        if qpoints:
            parts = qpoints.split()
            if len(parts) != 3 or any(int(x) <= 0 for x in parts):
                from core.exceptions import InvalidParameterError
                raise InvalidParameterError('qpoints', qpoints, '3个正整数')

        return True
