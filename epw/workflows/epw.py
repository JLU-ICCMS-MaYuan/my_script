"""
EPW电声耦合计算工作流

从epw_run.py重构而来
提供Electron-Phonon Wannier (EPW)计算功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import os
import logging

from epw.config import config
from epw.epw_inputpara import epw_inputpara
from epw.epw_writeinput import epw_writeinput
from epw.epw_writesubmit import epw_writesubmit
from epw.epw_submitjob import epw_submitjob

logger = logging.getLogger(__name__)


def check_pid_jobid(ids: list, submit_job_system):
    """检查任务是否还在运行"""
    if submit_job_system == "bash":
        i = 0
        while True:
            osawk = """sleep 5 | ps -ef | grep -E "pw.x" |  grep -v grep | awk '{print $2}'"""
            _jobids = os.popen(osawk).read()
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return
    elif submit_job_system == "slurm":
        i = 0
        while True:
            osawk = """sleep 5 | squeue | awk '{print $1}'"""
            _jobids = os.popen(osawk).read()
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return
    elif submit_job_system == "pbs":
        i = 0
        while True:
            osawk = """sleep 5 | qstat | awk '{print $1}' | cut -d . -f1"""
            _jobids = os.popen(osawk).read()
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return


class EPWWorkflow:
    """
    EPW电声耦合计算工作流

    支持8种计算模式：
    - epw_eband: Wannier拟合能带
    - epw_phono: 电声矩阵元傅里叶变换 + 插值
    - epw_phonodata: 仅插值电子/声子
    - epw_elph: 电声耦合计算
    - epw_sc: 超导计算（各向同性/各向异性）
    - epw_prtgkk: 输出电声矩阵元
    - epw_fermi_nest: 费米嵌套
    - epw_linearized_iso: 线性化各向同性计算
    """

    def __init__(self, args):
        # 读取配置
        _config = config(args).read_config()

        # 准备输入参数
        self.epw_inputpara = epw_inputpara.init_from_config(_config)
        self.epw_writeinput = epw_writeinput(self.epw_inputpara)
        self.epw_writesubmit = epw_writesubmit(self.epw_inputpara)
        self.epw_submitjob = epw_submitjob(self.epw_inputpara)

        self.run()

    def run(self):
        """根据模式执行相应的计算"""
        if self.epw_inputpara.mode == "epw_eband":
            logger.info("Perform an Wannier calculation and to fit energy bands")
            self.epw_eband()
        elif self.epw_inputpara.mode == "epw_phono":
            logger.info("Perform an EPW calculation to Fourier-transform the electron-phonon matrix element from a coarse k and q-point grids to real space and then interpolate the electronic band structure and phononic dispersion along the high symmetry line by reading modified_`prefix`_band.kpt.")
            self.epw_phono()
        elif self.epw_inputpara.mode == "epw_phonodata":
            logger.info("Only perform an EPW calculation to interpolate the electronic band structure and phononic dispersion along the high symmetry line by reading modified_`prefix`_band.kpt.")
            self.epw_phonodata()
        elif self.epw_inputpara.mode == "epw_elph":
            self.epw_elph()
        elif self.epw_inputpara.mode == "epw_sc":
            self.epw_sc()
        elif self.epw_inputpara.mode == "epw_prtgkk":
            self.epw_prtgkk()
        elif self.epw_inputpara.mode == "epw_fermi_nest":
            self.epw_fermi_nest()
        elif self.epw_inputpara.mode == "epw_linearized_iso":
            self.epw_linearized_iso()
        else:
            raise ValueError("Invalid mode selected.")

    def epw_eband(self):
        """Wannier拟合能带"""
        # 生成输入文件
        inputfilename = self.epw_writeinput.writeinput(mode="epw_eband")
        logger.info(inputfilename)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_eband")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename, jobname)

    def epw_phono(self):
        """电声矩阵元傅里叶变换 + 插值"""
        # 生成输入文件
        inputfilename1, inputfilename2 = self.epw_writeinput.writeinput(mode="epw_phono")
        logger.info(inputfilename1)
        logger.info(inputfilename2)
        print(inputfilename1, inputfilename2)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts([inputfilename1, inputfilename2], mode="epw_phono")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename1, jobname)

    def epw_phonodata(self):
        """仅插值电子/声子"""
        # 生成输入文件
        inputfilename1, inputfilename2 = self.epw_writeinput.writeinput(mode="epw_phonodata")
        logger.info(inputfilename1)
        logger.info(inputfilename2)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts([inputfilename1, inputfilename2], mode="epw_phonodata")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename1, jobname)

    def epw_elph(self):
        """电声耦合计算"""
        # 生成输入文件
        inputfilename = self.epw_writeinput.writeinput(mode="epw_elph")
        logger.info(inputfilename)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_elph")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename, jobname)

    def epw_sc(self):
        """超导计算（各向同性/各向异性）"""
        # 生成输入文件
        inputfilename1, inputfilename2 = self.epw_writeinput.writeinput(mode="epw_sc")
        logger.info(inputfilename1)
        logger.info(inputfilename2)

        # 生成提交脚本
        jobname1, jobname2 = self.epw_writesubmit.write_submit_scripts([inputfilename1, inputfilename2], mode="epw_sc")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode3([jobname1, jobname2])

    def epw_prtgkk(self):
        """输出电声矩阵元"""
        # 生成输入文件
        inputfilename = self.epw_writeinput.writeinput(mode="epw_prtgkk")
        logger.info(inputfilename)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_prtgkk")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode2(inputfilename, jobname, "prtgkk")

    def epw_fermi_nest(self):
        """费米嵌套计算"""
        # 生成输入文件
        inputfilename = self.epw_writeinput.writeinput(mode="epw_fermi_nest")
        logger.info(inputfilename)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_fermi_nest")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode2(inputfilename, jobname, "fermi_nest")

    def epw_linearized_iso(self):
        """线性化各向同性计算"""
        # 生成输入文件
        inputfilename1 = self.epw_writeinput.writeinput(mode="epw_linearized_iso")
        logger.info(inputfilename1)

        # 生成提交脚本
        jobname = self.epw_writesubmit.write_submit_scripts([inputfilename1, None], mode="epw_linearized_iso")

        # 提交任务
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode3([jobname, None])
