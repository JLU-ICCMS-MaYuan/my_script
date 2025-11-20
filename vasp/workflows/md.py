"""
VASP分子动力学工作流

从vasp_run.py重构而来
提供分子动力学模拟功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import logging
from argparse import ArgumentParser
from pathlib import Path

from vasp.config import config
from vasp.vasp_inputpara import vasp_mdpara
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class MDWorkflow:
    """
    VASP分子动力学工作流

    执行分子动力学模拟计算，包括：
    - 准备POSCAR和POTCAR
    - 生成MD相关的INCAR
    - 生成KPOINTS（均匀网格或Gamma点）
    - 提交任务
    """

    def __init__(self, args: ArgumentParser):
        # 读取配置
        _config = config(args).read_config()

        # 准备POSCAR POTCAR
        md_inputpara = vasp_mdpara.init_from_config1(_config)

        # 初始化INCAR
        self._vasp_writeincar = vasp_writeincar(md_inputpara)
        self._vasp_writeincar.writeinput()

        # 初始化KPOINTS
        if md_inputpara.kspacing is not None:
            supercell_lattice = md_inputpara.struct_type.lattice.matrix
            md_inputpara.write_evenly_kpoints(
                lattice=supercell_lattice,
                kspacing=md_inputpara.kspacing,
                kpoints_path=md_inputpara.work_path,
            )
        else:
            md_inputpara.write_gamma_kpoints(
                kpoints_path=md_inputpara.work_path,
            )

        # 初始化提交脚本
        _vasp_writesubmit = vasp_writesubmit(md_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()

        # 提交任务
        _vasp_submitjob = vasp_submitjob(md_inputpara)
        if md_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
