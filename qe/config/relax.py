"""
结构优化配置模块

作者：Claude
创建时间：2025-11-19
"""

from config.base import BaseConfig, DEFAULT_PARAMS


class RelaxConfig(BaseConfig):
    """结构优化配置"""

    def _set_defaults(self):
        """设置默认参数"""
        # 基础参数
        self.set_param('ecutwfc', DEFAULT_PARAMS['ecutwfc'])
        self.set_param('ecutrho', DEFAULT_PARAMS['ecutrho'])

        # 收敛参数
        self.set_param('conv_thr', '1.0d-8')
        self.set_param('forc_conv_thr', '1.0d-6')
        self.set_param('etot_conv_thr', '1.0d-7')

        # 计算模式
        self.set_param('calculation', 'relax')  # or 'vc-relax'
        self.set_param('nstep', 100)

        # K点
        self.set_param('kpoints', '6 6 6')

    def validate(self) -> bool:
        """验证参数"""
        # 验证ecutwfc范围
        ecutwfc = self.get_param('ecutwfc')
        if ecutwfc <= 0 or ecutwfc > 200:
            from core.exceptions import InvalidParameterError
            raise InvalidParameterError('ecutwfc', ecutwfc, '0 < ecutwfc <= 200')

        return True
