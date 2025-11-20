"""
VASP Readers模块

整合了原vasptools的功能，提供统一的VASP输出文件读取接口：
- poscar: POSCAR文件读取
- outcar: OUTCAR文件读取
- vasprunxml: vasprun.xml文件读取
- parser_vasp: VASP输出解析

作者：Claude (整合自vasptools)
创建时间：2025-11-20
"""

# 保持向后兼容，重新导出vasptools中的所有功能
try:
    from .vasptools.poscar import *
    from .vasptools.outcar import *
    from .vasptools.vasprunxml import *
    from .vasptools.parser_vasp import *
except ImportError:
    # 如果导入失败，提供一个友好的错误提示
    import warnings
    warnings.warn("vasptools模块未完全加载，某些功能可能不可用")

__all__ = []
