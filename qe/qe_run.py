import time
import os
import logging
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain

from config import config
from qe_inputpara import * 
from qe_writeinput import qe_writeinput
from qe_writesubmit import qe_writesubmit
from qe_submitjob import qe_submitjob


def check_pid_jobid(ids: list, submit_job_system):
    if submit_job_system == "bash":
        i = 0
        while True:
            osawk = """sleep 5 | ps -ef | grep -E "pw.x" |  grep -v grep | awk '{print $2}'""" 
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
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
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return



class qe_relax:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.relax_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the input
        self.qe_writeinput  = qe_writeinput.init_from_relaxinput(self.relax_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_relaxinput(self.relax_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
        
        # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
        self.qe_submitjob = qe_submitjob.init_from_relaxinput(self.relax_inputpara)
        if self.relax_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_scf:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.scf_inputpara  = qe_inputpara.init_from_config(self._config)

        # init the input
        self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.scf_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_scfinput(self.scf_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

        # submit the job
        self.qe_submitjob = qe_submitjob.init_from_scfinput(self.scf_inputpara)
        if self.scf_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_phono:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.phono_inputpara = qephono_inputpara.init_from_config(self._config)

        if self.phono_inputpara.mode == "merge":
            self.phono_inputpara.merge(self.phono_inputpara.work_underpressure)
        else:
            # init the input
            self.qe_writeinput  = qe_writeinput.init_from_phonoinput(self.phono_inputpara)
            inputfilename = self.qe_writeinput.writeinput()

            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_phonoinput(self.phono_inputpara)
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename)

            # submit the job
            self.qe_submitjob = qe_submitjob.init_from_phonoinput(self.phono_inputpara)
            if self.phono_inputpara.queue is not None :
                if self.phono_inputpara.mode == "nosplit":
                    self.qe_submitjob.submit_mode2(inputfilename, jobnames)
                elif self.phono_inputpara.mode == "split_dyn0":
                    self.qe_submitjob.submit_mode3(inputfilename, jobnames)
                elif self.phono_inputpara.mode == "split_assignQ":
                    self.qe_submitjob.submit_mode4(inputfilename, jobnames)
                else:
                    self.qe_submitjob.submit_mode1(inputfilename, jobnames)
        

class qe_dos:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.dos_inputpara = qedos_inputpara.init_from_config(self._config)

        # write input parameter
        self.qe_writeinput  = qe_writeinput.init_from_dosinput(self.dos_inputpara)
        inputfilename = self.qe_writeinput.writeinput()
        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_dosinput(self.dos_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)
        # submit the job
        self.qe_submitjob = qe_submitjob.init_from_dosinput(self.dos_inputpara)
        if self.dos_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_superconduct:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare input parameter
        self.sc_inputpara  = qesc_inputpara.init_from_config(self._config)
        # init the input

        self.qe_writeinput  = qe_writeinput.init_from_scinput(self.sc_inputpara)
        inputfilename = self.qe_writeinput.writeinput()

        # init the submit job script
        self.qe_writesubmit = qe_writesubmit.init_from_scinput(self.sc_inputpara)
        jobname = self.qe_writesubmit.write_submit_scripts(inputfilename)

        # submit the job
        self.qe_submitjob   = qe_submitjob.init_from_scinput(self.sc_inputpara)
        if self.sc_inputpara.queue is not None:
            self.qe_submitjob.submit_mode1(inputfilename, jobname)


class qe_prepare:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        #进行结构弛豫
        self.prepare_inputpara  = qeprepare_inputpara.init_from_config(self._config)
        if "relax-vc" in self.prepare_inputpara.mode:
            self.qe_writeinput  = qe_writeinput.init_from_relaxinput(self.prepare_inputpara)
            inputfilename = self.qe_writeinput.writeinput(mode="relax-vc")
            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_relaxinput(self.prepare_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename, mode="relax-vc")
            
            # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
            self.qe_submitjob = qe_submitjob.init_from_relaxinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None:
                ids = self.qe_submitjob.submit_mode1(inputfilename, jobname)
            logger.info("The program is running relax")
            try:
                check_pid_jobid(ids, self.qe_submitjob.submit_job_system)
            except:
                logger.info(f"ids={None}")
            logger.info("The program finished relax")
            ids = []

        # 读入结构弛豫的输出文件进行自洽计算和声子计算
        if self.prepare_inputpara.input_file_path.name != "relax.out":
            input_file_path = Path(self.prepare_inputpara.work_underpressure).joinpath("relax.out")
        else:
            input_file_path = self.prepare_inputpara.input_file_path
        if not input_file_path.exists():
            raise FileExistsError(f"{input_file_path.absolute()} doesn't exist!!! ")
        self.prepare_inputpara  = qeprepare_inputpara(
            work_path=self.prepare_inputpara.work_path,
            press=self.prepare_inputpara.press,
            submit_job_system=self.prepare_inputpara.submit_job_system,
            input_file_path=input_file_path,
            **self._config,)
        if "scffit" in self.prepare_inputpara.mode:
            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename = self.qe_writeinput.writeinput(mode="scffit")
            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_scfinput(self.prepare_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename, mode="scffit")
            
            # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
            self.qe_submitjob = qe_submitjob.init_from_scfinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None:
                ids = self.qe_submitjob.submit_mode1(inputfilename, jobname)
            logger.info("The program is running scffit")
            try:
                check_pid_jobid(ids, self.qe_submitjob.submit_job_system)
            except:
                logger.info(f"ids={None}")
            logger.info("The program finished scffit")
            ids = []

        if "scf" in self.prepare_inputpara.mode:
            self.qe_writeinput  = qe_writeinput.init_from_scfinput(self.prepare_inputpara)
            inputfilename = self.qe_writeinput.writeinput(mode="scf")
            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_scfinput(self.prepare_inputpara)
            jobname = self.qe_writesubmit.write_submit_scripts(inputfilename, mode="scf")
            
            # submit the job. If we didn't set the parameter of `queue`, it will be set `None` in `qe_inputpara`
            self.qe_submitjob = qe_submitjob.init_from_scfinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None:
                ids = self.qe_submitjob.submit_mode1(inputfilename, jobname)
            logger.info("The program is running scf")
            try:
                check_pid_jobid(ids, self.qe_submitjob.submit_job_system)
            except:
                logger.info(f"ids={None}")
            logger.info("The program finished scf")
            ids = []
        
        if "nosplit" in self.prepare_inputpara.mode and self.prepare_inputpara.dyn0_flag:
            self.qe_writeinput  = qe_writeinput.init_from_phonoinput(self.prepare_inputpara)
            inputfilename = self.qe_writeinput.writeinput(mode="nosplit")

            # init the submit job script
            self.qe_writesubmit = qe_writesubmit.init_from_phonoinput(self.prepare_inputpara)
            jobnames = self.qe_writesubmit.write_submit_scripts(inputfilename, mode="nosplit")

            # submit the job
            self.qe_submitjob = qe_submitjob.init_from_phonoinput(self.prepare_inputpara)
            if self.prepare_inputpara.queue is not None :
                self.qe_submitjob.submit_mode2(inputfilename, jobnames)
            

