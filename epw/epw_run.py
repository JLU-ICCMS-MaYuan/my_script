import os
import logging
from argparse import ArgumentParser

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


def epw_band(args):
    # read input para
    _config = config(args).read_config()

    # prepare input parameter
    epw_inputpara = epw_inputpara.init_from_config(_config)

    #  # init the input
    epw_writeinput = epw_writeinput(epw_inputpara)
    inputfilename = epw_writeinput.writeinput()
    logger.info(inputfilename)
    
    # init the submit job script
    epw_writesubmit = epw_writesubmit(epw_inputpara)
    jobname = epw_writesubmit.write_submit_scripts(inputfilename)

    # submit the job
    epw_submitjob = epw_submitjob(epw_inputpara)
    if epw_inputpara.queue is not None:
        epw_submitjob.submit_mode1(inputfilename, jobname)