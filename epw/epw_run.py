import os
import logging

from epw.config import config
from epw.epw_inputpara import epw_inputpara
from epw.epw_writeinput import epw_writeinput
from epw.epw_writesubmit import epw_writesubmit
from epw.epw_submitjob import epw_submitjob 

logger = logging.getLogger(__name__)

logger.info("This is an info message")
logger.debug("This is a debug message")


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
    elif submit_job_system == "pbs":
        i = 0
        while True:
            osawk = """sleep 5 | qstat | awk '{print $1}' | cut -d . -f1""" 
            _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
            jobids = _jobids.strip("\n").split("\n")
            for id in ids:
                if id not in jobids:
                    i += 1
                if i == len(ids):
                    return


class epw_run:
    def __init__(self, args):
        # read input para
        _config = config(args).read_config()

        # prepare input parameter
        self.epw_inputpara = epw_inputpara.init_from_config(_config)
        self.epw_writeinput = epw_writeinput(self.epw_inputpara)
        self.epw_writesubmit = epw_writesubmit(self.epw_inputpara)
        self.epw_submitjob = epw_submitjob(self.epw_inputpara)
        
        self.run()
        
    def run(self):
        if self.epw_inputpara.mode == "epw_eband":
            self.epw_eband()
        elif self.epw_inputpara.mode == "epw_phono":
            self.epw_phono()
        elif self.epw_inputpara.mode == "epw_elph":
            self.epw_elph()
        else:
            raise ValueError("Invalid mode selected.")
    def epw_eband(self, _mode):

        # init the input file
        inputfilename = self.epw_writeinput.writeinput(mode="epw_eband")
        logger.info(inputfilename)
        
        # init the submit job script
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_eband")

        # submit the job
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename, jobname)
    
    def epw_phono(self):
        # init the input file
        inputfilename1, inputfilename2, inputfilename3 = self.epw_writeinput.writeinput(mode="epw_phono")
        logger.info(inputfilename1)
        logger.info(inputfilename2)
        logger.info(inputfilename3)
        print(inputfilename1, inputfilename2, inputfilename3)
        # init the submit job script
        jobname = self.epw_writesubmit.write_submit_scripts([inputfilename1, inputfilename2, inputfilename3], mode="epw_phono")

        # submit the job
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename1, jobname)
            
    def epw_elph(self):
        # init the input file
        inputfilename = self.epw_writeinput.writeinput(mode="epw_elph")
        logger.info(inputfilename)
        
        # init the submit job script
        jobname = self.epw_writesubmit.write_submit_scripts(inputfilename, mode="epw_elph")
        
        # submit the job    
        if self.epw_inputpara.queue is not None:
            self.epw_submitjob.submit_mode1(inputfilename, jobname)