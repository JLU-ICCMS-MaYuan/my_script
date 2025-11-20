"""
VASP结构弛豫工作流

从vasp_run.py重构而来
提供结构优化计算功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import logging
from argparse import ArgumentParser
from pathlib import Path

from vasp.config import config
from vasp.vasp_inputpara import vasp_inputpara
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class RelaxWorkflow:
    """
    VASP结构弛豫工作流

    执行结构优化计算，包括：
    - 准备POSCAR和POTCAR
    - 生成INCAR输入文件
    - 编写作业提交脚本
    - 提交任务到队列系统
    """

    def __init__(self, args: ArgumentParser):
        # 读取配置
        _config = config(args).read_config()

        # 准备POSCAR POTCAR
        self.relax_inputpara = vasp_inputpara.init_from_config1(_config)

        # 初始化INCAR
        self._vasp_writeincar = vasp_writeincar(self.relax_inputpara)
        self._vasp_writeincar.writeinput()

        # 初始化提交脚本
        _vasp_writesubmit = vasp_writesubmit(self.relax_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()

        # 提交任务
        _vasp_submitjob = vasp_submitjob(self.relax_inputpara)
        if self.relax_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
