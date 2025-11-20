"""
Common Core模块

提供QE、VASP、EPW三个项目共享的核心功能：
- 配置读取 (config.py)
- 任务提交 (submit.py)
- 日志系统 (logging.py)
- 基础类 (base.py)

作者：Claude
创建时间：2025-11-20
"""

from .config import (
    BaseConfig,
    QEConfig,
    VASPConfig,
    EPWConfig,
    create_config,
)

from .submit import (
    QueueSystem,
    JobSubmitter,
    submit_job_simple,
)

from .logging import (
    setup_logger,
    setup_default_logger,
    LoggerAdapter,
    get_qe_logger,
    get_vasp_logger,
    get_epw_logger,
)

from .base import (
    BaseWorkflow,
    BaseInputWriter,
    BaseOutputReader,
    StructureFile,
)

__all__ = [
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
    'LoggerAdapter',
    'get_qe_logger',
    'get_vasp_logger',
    'get_epw_logger',

    # Base
    'BaseWorkflow',
    'BaseInputWriter',
    'BaseOutputReader',
    'StructureFile',
]
