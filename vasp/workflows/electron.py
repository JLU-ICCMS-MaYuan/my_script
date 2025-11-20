"""
VASP电子结构计算工作流

从vasp_run.py重构而来
提供电子能带、DOS、COHP等计算功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import os
import sys
import shutil
import logging
from argparse import ArgumentParser
from pathlib import Path

from vasp.config import config
from vasp.vasp_inputpara import vasp_eletronpara
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class ElectronWorkflow:
    """
    VASP电子结构计算工作流

    支持多种计算模式：
    - scf: 自洽场计算
    - eband: 电子能带
    - eledos: 电子态密度
    - cohp: 晶体轨道哈密顿布居

    可以组合计算（如scf+eband+eledos）或单独计算
    """

    def __init__(self, args: ArgumentParser):
        # 读取配置
        _config = config(args).read_config()

        # 准备POSCAR POTCAR
        self.eletron_inputpara = vasp_eletronpara.init_from_config1(_config)

        # TODO 简化整个计算流程
        # 准备输入文件
        self._vasp_writeincar = vasp_writeincar(self.eletron_inputpara)

        if 'scf' in self.eletron_inputpara.mode:
            scf_path = self.scf(self.eletron_inputpara.kspacing)
        if 'eband' in self.eletron_inputpara.mode:
            eband_path = self.eband()
        if 'eledos' in self.eletron_inputpara.mode:
            eledos_path = self.eledos(self.eletron_inputpara.kspacing / 2)
        if 'cohp' in self.eletron_inputpara.mode:
            cohp_path = self.cohp(self.eletron_inputpara.kspacing)

        # 准备提交任务的脚本
        self._prepare_and_submit_jobs(
            scf_path if 'scf' in self.eletron_inputpara.mode else None,
            eband_path if 'eband' in self.eletron_inputpara.mode else None,
            eledos_path if 'eledos' in self.eletron_inputpara.mode else None,
            cohp_path if 'cohp' in self.eletron_inputpara.mode else None
        )

    def _prepare_and_submit_jobs(self, scf_path, eband_path, eledos_path, cohp_path):
        """准备和提交各种组合的计算任务"""
        mode = self.eletron_inputpara.mode

        # 同时进行计算的任务, 作业脚本放在work_path
        if 'scf' in mode and 'eledos' in mode and 'eband' in mode:
            self._submit_combined_job("scf-eband-eledos")
        elif 'scf' in mode and 'eband' in mode and 'eledos' not in mode:
            self._submit_combined_job("scf-eband")
        elif 'scf' in mode and 'eledos' in mode and 'eband' not in mode:
            self._submit_combined_job("scf-eledos")
        elif 'scf' not in mode and 'eledos' in mode and 'eband' in mode:
            self._submit_combined_job("eband-eledos", need_chgcar=True)

        # 单独进行计算的任务, 作业脚本放在work_path/sub_workpath
        elif 'scf' in mode and 'eledos' not in mode and 'eband' not in mode:
            self._submit_single_job("scf", scf_path)
        elif 'scf' not in mode and 'eledos' not in mode and 'eband' in mode:
            self._submit_single_job("eband", eband_path, need_chgcar=True)
        elif 'scf' not in mode and 'eledos' in mode and 'eband' not in mode:
            self._submit_single_job("eledos", eledos_path, need_chgcar=True)

        # cohp计算
        if 'cohp' in mode:
            self._submit_single_job("cohp", cohp_path)

    def _submit_combined_job(self, mode, need_chgcar=False):
        """提交组合任务（在work_path提交）"""
        if need_chgcar:
            chgcar_src = self.eletron_inputpara.work_path.joinpath("scf", "CHGCAR")
            if not chgcar_src.exists():
                logger.error(f"The CHGCAR is not found in path \n{chgcar_src.absolute()}")
                logger.error("So The program will exit")
                sys.exit(1)

        _vasp_writesubmit = vasp_writesubmit(self.eletron_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts(mode=mode)

        _vasp_submitjob = vasp_submitjob(self.eletron_inputpara)
        if self.eletron_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)

    def _submit_single_job(self, mode, submit_path, need_chgcar=False):
        """提交单独任务（在子目录提交）"""
        if need_chgcar:
            chgcar_src = self.eletron_inputpara.work_path.joinpath("scf", "CHGCAR")
            chgcar_dst = submit_path.joinpath("CHGCAR")
            if not chgcar_src.exists():
                logger.error(f"The CHGCAR is not found in path \n{chgcar_src.absolute()}")
                logger.error("So The program will exit")
                sys.exit(1)
            else:
                shutil.copy(chgcar_src, chgcar_dst)

        _vasp_writesubmit = vasp_writesubmit(self.eletron_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts(mode=mode, submitjob_path=submit_path)

        _vasp_submitjob = vasp_submitjob(self.eletron_inputpara)
        if self.eletron_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname, submit_path=submit_path)

    def scf(self, kspacing):
        """准备SCF自洽场计算"""
        scf_path = self.eletron_inputpara.work_path.joinpath("scf")
        if not scf_path.exists():
            os.mkdir(scf_path)

        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            scf_path,
        )
        self.eletron_inputpara.get_potcar(scf_path)

        # 准备计算电子自洽的INCAR
        self._vasp_writeincar.writeinput(mode="scf", incar_path=scf_path)

        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing,
            scf_path,
        )
        return scf_path

    def eband(self):
        """准备电子能带计算"""
        scf_path = self.eletron_inputpara.work_path.joinpath("scf")
        eband_path = self.eletron_inputpara.work_path.joinpath('eband')
        if not scf_path.exists():
            os.mkdir(scf_path)
        if not eband_path.exists():
            os.mkdir(eband_path)

        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            eband_path)

        # 准备POTCAR
        self.eletron_inputpara.get_potcar(eband_path)

        # 为计算电子band准备高对称路径
        self.eletron_inputpara.write_highsymmetry_kpoints(
            self.eletron_inputpara.ase_type,
            kpoints_path=eband_path,
            autoselect=self.eletron_inputpara.autoselect,
            vaspkitflag=self.eletron_inputpara.vaspkitflag,
        )

        # 准备计算电子band的INCAR
        self._vasp_writeincar.writeinput(mode="eband", incar_path=eband_path)
        return eband_path

    def eledos(self, kspacing):
        """准备电子态密度计算"""
        logger.debug("KSPACING in `eledos` have to be twice than that in `scf` ")
        scf_path = self.eletron_inputpara.work_path.joinpath("scf")
        eledos_path = self.eletron_inputpara.work_path.joinpath('eledos')
        if not scf_path.exists():
            os.mkdir(scf_path)
        if not eledos_path.exists():
            os.mkdir(eledos_path)

        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            eledos_path)
        self.eletron_inputpara.get_potcar(eledos_path)

        # 准备计算电子DOS的INCAR
        self._vasp_writeincar.writeinput(mode='eledos', incar_path=eledos_path)

        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing,
            kpoints_path=eledos_path,
        )
        return eledos_path

    def cohp(self, kspacing):
        """准备COHP计算"""
        cohp_path = self.eletron_inputpara.work_path.joinpath('cohp')
        if not cohp_path.exists():
            os.mkdir(cohp_path)

        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            cohp_path)
        self.eletron_inputpara.get_potcar(cohp_path)

        # 准备计算电子cohp的INCAR
        self._vasp_writeincar.writeinput(mode='cohp', incar_path=cohp_path)

        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing,
            kpoints_path=cohp_path,
        )
        return cohp_path
