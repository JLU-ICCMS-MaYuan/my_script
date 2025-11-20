"""
VASP IO模块

提供VASP文件读取和写入功能：
- readers: 读取VASP输出文件（整合了vasptools）
- writers: 写入VASP输入文件

作者：Claude
创建时间：2025-11-20
"""

from . import readers
from . import writers

__all__ = ['readers', 'writers']
