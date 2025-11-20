"""
VASP声子计算工作流

从vasp_run.py重构而来
提供声子谱和声子DOS计算功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import logging
from argparse import ArgumentParser
from pathlib import Path

from vasp.config import config
from vasp.vasp_inputpara import vasp_phonopara
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class PhononWorkflow:
    """
    VASP声子计算工作流

    支持两种计算模式：
    - disp: 位移超胞方法
    - dfpt: 密度泛函微扰理论

    包括：
    - 准备POSCAR和POTCAR
    - 生成INCAR和KPOINTS
    - 提交任务
    """

    def __init__(self, args: ArgumentParser):
        # 读取配置
        _config = config(args).read_config()

        # 准备POSCAR POTCAR
        phono_inputpara = vasp_phonopara.init_from_config1(_config)

        # 初始化INCAR
        self._vasp_writeincar = vasp_writeincar(phono_inputpara)
        self._vasp_writeincar.writeinput()

        # 初始化KPOINTS
        if phono_inputpara.kdensity is not None:
            phono_inputpara.create_kpoints_by_pymatgen(
                phono_inputpara.sposcar_struct_type,
                phono_inputpara.work_path.joinpath("KPOINTS"),
                phono_inputpara.kdensity,
            )
        elif phono_inputpara.kspacing is not None:
            supercell_lattice = phono_inputpara.sposcar_struct_type.lattice.matrix
            phono_inputpara.write_evenly_kpoints(
                lattice=supercell_lattice,
                kspacing=phono_inputpara.kspacing,
                kpoints_path=phono_inputpara.work_path,
            )

        # 初始化提交脚本
        _vasp_writesubmit = vasp_writesubmit(phono_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()

        # 提交任务
        _vasp_submitjob = vasp_submitjob(phono_inputpara)
        if phono_inputpara.queue is not None:
            if phono_inputpara.mode == 'disp':
                _vasp_submitjob.submit_mode2(jobname)
            elif phono_inputpara.mode == 'dfpt':
                _vasp_submitjob.submit_mode1(jobname)
