"""
Common - 计算软件共享核心库

为QE、VASP、EPW三个计算软件提供统一的核心功能。
消除80%+的代码重复，提供一致的接口和行为。

主要模块：
- core: 核心功能（配置、提交、日志、基类）
- utils: 工具函数（结构处理、文件操作）

作者：Claude
创建时间：2025-11-20
版本：1.0.0
"""

__version__ = '1.0.0'

# 导出核心模块
from . import core
from . import utils

# 便捷导入
from .core import (
    # Config
    BaseConfig,
    QEConfig,
    VASPConfig,
    EPWConfig,
    create_config,

    # Submit
    QueueSystem,
    JobSubmitter,
    submit_job_simple,

    # Logging
    setup_logger,
    setup_default_logger,
    get_qe_logger,
    get_vasp_logger,
    get_epw_logger,

    # Base
    BaseWorkflow,
    BaseInputWriter,
    BaseOutputReader,
    StructureFile,
)

from .utils import (
    # Structure
    read_poscar,
    write_poscar,
    poscar_to_qe_format,

    # File Ops
    ensure_dir,
    copy_file,
    remove_file,
    clean_tmp_files,
    find_files,
)

__all__ = [
    'core',
    'utils',

    # Config
    'BaseConfig',
    'QEConfig',
    'VASPConfig',
    'EPWConfig',
    'create_config',

    # Submit
    'QueueSystem',
    'JobSubmitter',
    'submit_job_simple',

    # Logging
    'setup_logger',
    'setup_default_logger',
    'get_qe_logger',
    'get_vasp_logger',
    'get_epw_logger',

    # Base
    'BaseWorkflow',
    'BaseInputWriter',
    'BaseOutputReader',
    'StructureFile',

    # Structure
    'read_poscar',
    'write_poscar',
    'poscar_to_qe_format',

    # File Ops
    'ensure_dir',
    'copy_file',
    'remove_file',
    'clean_tmp_files',
    'find_files',
]
