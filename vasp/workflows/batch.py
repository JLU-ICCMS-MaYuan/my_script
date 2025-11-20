"""
VASP批量计算工作流

从vasp_run.py重构而来
提供批量结构计算和批量声子计算功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import os
import sys
import logging
from argparse import ArgumentParser
from pathlib import Path

from vasp.config import config
from vasp.vasp_inputpara import vaspbatch_inputpara, vaspbatch_phonopara
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class BatchWorkflow:
    """
    VASP批量计算工作流

    批量处理多个结构文件（.cif, .vasp, .res）
    支持多个压力点计算
    """

    def __init__(self, args: ArgumentParser):
        _config = config(args).read_config()
        input_dir_path = Path(_config['input_file_path'])

        if not input_dir_path.is_dir():
            logger.error(f"The {input_dir_path} path doesn't exist!")
            sys.exit(1)

        # 获取所有结构文件
        input_files_path = [
            input_dir_path.joinpath(pathname)
            for pathname in os.listdir(input_dir_path)
            if input_dir_path.joinpath(pathname).suffix in [".cif", ".vasp", ".res"]
        ]

        if not input_files_path:
            logger.error(f"The program didn't get any structures from {input_dir_path}")
            sys.exit(1)

        # 提取配置参数
        work_path = _config['work_path']
        del _config['work_path']
        presses = _config['presses']
        del _config['presses']
        press = _config['press']
        del _config['press']
        submit_job_system = _config['submit_job_system']
        del _config['submit_job_system']
        del _config['input_file_path']
        pp_dir = _config['pp_dir']
        del _config['pp_dir']

        # 批量处理每个结构文件
        for input_file_path in input_files_path:
            logger.info(f"Create directory for {input_file_path} file !!!")

            if presses is not None:
                # 多个压力点
                for press in presses:
                    logger.debug(press)
                    self._process_single_structure(
                        work_path, press, submit_job_system,
                        input_file_path, pp_dir, _config
                    )
            else:
                # 单个压力点
                self._process_single_structure(
                    work_path, None, submit_job_system,
                    input_file_path, pp_dir, _config
                )

    def _process_single_structure(self, work_path, press, submit_job_system,
                                   input_file_path, pp_dir, config_dict):
        """处理单个结构"""
        self.batch_inputpara = vaspbatch_inputpara(
            work_path=work_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            pp_dir=pp_dir,
            **config_dict
        )

        # 初始化INCAR
        self._vasp_writeincar = vasp_writeincar(self.batch_inputpara)
        self._vasp_writeincar.writeinput()

        # 初始化提交脚本
        _vasp_writesubmit = vasp_writesubmit(self.batch_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()

        # 提交任务
        _vasp_submitjob = vasp_submitjob(self.batch_inputpara)
        if self.batch_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)


class BatchPhononWorkflow:
    """
    VASP批量声子计算工作流

    批量处理多个.vasp结构文件
    执行声子计算
    """

    def __init__(self, args: ArgumentParser):
        _config = config(args).read_config()
        input_dir_path = Path(_config['input_file_path'])

        if not input_dir_path.is_dir():
            logger.error(f"The {input_dir_path} path doesn't exist!")
            sys.exit(1)

        input_files_path = list(input_dir_path.glob("*.vasp"))

        # 提取配置参数
        work_path = _config['work_path']
        del _config['work_path']
        press = _config['press']
        del _config['press']
        submit_job_system = _config['submit_job_system']
        del _config['submit_job_system']
        del _config['input_file_path']
        pp_dir = _config['pp_dir']
        del _config['pp_dir']

        for input_file_path in input_files_path:
            # 准备POSCAR POTCAR
            phono_inputpara = vaspbatch_phonopara(
                work_path=work_path,
                press=press,
                submit_job_system=submit_job_system,
                input_file_path=input_file_path,
                pp_dir=pp_dir,
                **_config
            )

            # 初始化INCAR
            self._vasp_writeincar = vasp_writeincar(phono_inputpara)
            self._vasp_writeincar.writeinput()

            # 初始化KPOINTS
            phono_inputpara.create_kpoints_by_pymatgen(
                phono_inputpara.sposcar_struct_type,
                phono_inputpara.work_path.joinpath("KPOINTS"),
                phono_inputpara.kdensity,
            )

            # 初始化提交脚本
            _vasp_writesubmit = vasp_writesubmit(phono_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts()

            # 提交任务
            _vasp_submitjob = vasp_submitjob(phono_inputpara)
            if phono_inputpara.queue is not None:
                _vasp_submitjob.submit_mode2(jobname)
