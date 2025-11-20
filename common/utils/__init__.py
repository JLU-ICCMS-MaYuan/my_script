"""
Common Utils模块

提供通用工具函数：
- 结构文件处理 (structure.py)
- 文件操作 (file_ops.py)

作者：Claude
创建时间：2025-11-20
"""

from .structure import (
    read_poscar,
    write_poscar,
    poscar_to_qe_format,
    get_lattice_parameters,
    direct_to_cartesian,
    cartesian_to_direct,
)

from .file_ops import (
    ensure_dir,
    copy_file,
    copy_files,
    move_file,
    remove_file,
    remove_dir,
    clean_tmp_files,
    find_files,
    get_latest_file,
    backup_file,
)

__all__ = [
    # Structure
    'read_poscar',
    'write_poscar',
    'poscar_to_qe_format',
    'get_lattice_parameters',
    'direct_to_cartesian',
    'cartesian_to_direct',

    # File Ops
    'ensure_dir',
    'copy_file',
    'copy_files',
    'move_file',
    'remove_file',
    'remove_dir',
    'clean_tmp_files',
    'find_files',
    'get_latest_file',
    'backup_file',
]
