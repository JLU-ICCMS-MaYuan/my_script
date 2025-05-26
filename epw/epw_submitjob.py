import re
import os
import sys
import time
import logging


from pathlib import Path

from epw.epw_inputpara import epw_inputpara
from epw.epwbin import epwbin_path

logger = logging.getLogger("epw_submitjob")

class epw_submitjob:
    
    def __init__(
        self,
        epw_inputpara: epw_inputpara
        ):

        self.epw_inputpara = epw_inputpara

        if self.epw_inputpara.submit_job_system == "slurm":
            self.submit_order = "sbatch"
        elif self.epw_inputpara.submit_job_system == "pbs":
            self.submit_order = "qsub"
        elif self.epw_inputpara.submit_job_system == "lsf":
            self.submit_order = "bsub <"
        elif self.epw_inputpara.submit_job_system == "bash":
            self.submit_order = "bash"
        else:
            self.submit_order = ''


    def submit_mode1(self, inputfilename, jobname):
        """
        submit_mode1 can be used to submit job by job-system, 
            such as by PBS, SLURM or Background executed mode of bash 
        """
        input_file = Path(self.epw_inputpara.work_path).joinpath(inputfilename)
        if not input_file.exists():
            raise FileExistsError(f" {inputfilename} doesn't exist")
        job_file = Path(self.epw_inputpara.work_path).joinpath(jobname)
        if not job_file.exists():
            raise FileExistsError(f" {jobname} doesn't exist")
        cwd = input_file.cwd()
        dst_dir = input_file.parent.absolute()
        os.chdir(dst_dir)
        if self.epw_inputpara.submit_job_system == "bash":
            logger.debug(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &")
            res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
            jobids = self.getpid()
        else:
            logger.debug(f"{self.submit_order} {jobname}")
            res = os.popen(f"{self.submit_order} {jobname}").read()
            jobids = re.findall(r"\d+", res)

        logger.info(f"{jobname} is executed. pid or jobids = {jobids}")
        os.chdir(cwd)
        # 检查任务是否成功提交，成功提交的话，应该会有进程号或者任务号返回。
        # 如果没有成功提交任务就跳出程序
        if not jobids:
            raise ValueError(f"The *.x of epw didn't run ! Because the jobids={jobids}. The program will exit! The order you use is {self.submit_order} {jobname}")
        return jobids
    
    def submit_mode2(self, inputfilename, jobname, submit_task_path):
        """专门用来提交计算prtgkk任务和fermi_nesting任务"""
        if not submit_task_path.exists():
            raise FileExistsError(f" {submit_task_path} doesn't exist")
        if not submit_task_path.joinpath(inputfilename).exists():
            raise FileExistsError(f" {inputfilename} doesn't exist")
        if not submit_task_path.joinpath(jobname).exists():
            raise FileExistsError(f" {jobname} doesn't exist")
        
        cwd = submit_task_path.cwd()
        dst_dir = submit_task_path.parent.absolute()
        os.chdir(dst_dir)
        if self.epw_inputpara.submit_job_system == "bash":
            logger.debug(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &")
            res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
            jobids = self.getpid()
        else:
            logger.debug(f"{self.submit_order} {jobname}")
            res = os.popen(f"{self.submit_order} {jobname}").read()
            jobids = re.findall(r"\d+", res)

        logger.info(f"{jobname} is executed. pid or jobids = {jobids}")
        os.chdir(cwd)
        # 检查任务是否成功提交，成功提交的话，应该会有进程号或者任务号返回。
        # 如果没有成功提交任务就跳出程序
        if not jobids:
            raise ValueError(f"The *.x of epw didn't run ! Because the jobids={jobids}. The program will exit! The order you use is {self.submit_order} {jobname}")
        return jobids

    def submit_mode3(self, inputfilename, jobname):
        """专门用来提交计算超导任务"""
        jobids = []
        iso_jobname = jobname[0]
        aniso_jobname = jobname[1]
        for mu in self.epw_inputpara.muc:
            iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(mu))
            aniso_mu_path = self.epw_inputpara.work_path.joinpath("aniso_muc_{}".format(mu))
            cwd = os.getcwd()
            # 提交iso的超导计算
            os.chdir(iso_mu_path)
            if self.epw_inputpara.submit_job_system == "bash":
                logger.debug(f"nohup {self.submit_order} {iso_jobname} > bash.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {iso_jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                logger.debug(f"{self.submit_order} {iso_jobname}")
                res = os.popen(f"{self.submit_order} {iso_jobname}").read()
                jobids = re.findall(r"\d+", res)
            logger.info(f"finish submit {iso_jobname}, jobids = {' '.join(jobids)}")
            os.chdir(cwd)
            
            # 提交aniso的超导计算
            os.chdir(aniso_mu_path) 
            if self.epw_inputpara.submit_job_system == "bash":
                logger.debug(f"nohup {self.submit_order} {aniso_jobname} > bash.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {aniso_jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                logger.debug(f"{self.submit_order} {aniso_jobname}")
                res = os.popen(f"{self.submit_order} {aniso_jobname}").read()
                jobids = re.findall(r"\d+", res)
            logger.info(f"finish submit {aniso_jobname}, jobids = {' '.join(jobids)}")
            os.chdir(cwd)

    @staticmethod
    def getpid():
        """get pid number"""
        jobids = []
        logger.debug("wait 6s, The program will tell you PID"); #time.sleep(6); 
        osawk = """ps -ef | grep -E "pw.x|ph.x|matdyn.x|lambda.x|q2r.x|eliashberg.x|dos.x|pp.x|projwfc.x" |  grep -v grep | awk '{print $2}'""" # return a series of number, such as: 423423 324233 423424
        # ps -ef ps -ef用于查看全格式的全部进程，其中“ps”是在Linux中是查看进程的命令，“-e ”参数代表显示所有进程，“-f”参数代表全格式。
        # grep -E  ‘grep’ ‘-E’ 选项表示使用扩展的正则表达式。如果你使用 ‘grep’ 命令时带 ‘-E’，你只需要用途 ‘|’ 来分隔OR条件。 grep -E 'pattern1|pattern2' filename
        # grep -v grep 这里可以比较看出，多出了一个进程号，这是grep时所多出来的进程，通过grep -v grep去除包含grep文本的进程行 ，避免影响最终数据的正确性
        #  awk '{print $2}' 这样，就可以抓取PID号
        _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
        jobids = _jobids.strip("\n").split("\n")
        return jobids

    @staticmethod
    def checksuffix(directory_path, suffix):
        """
        Checks for the existence of a file with a suffix in the destination directory
        directory_path : destination directory
        suffix : a file with a suffix

        return: True of False
        """ 

        files_underdir = os.listdir(directory_path)
        for file in files_underdir:
            if Path(file).suffix == suffix:
                return True

    @staticmethod
    def checkerror(directory_path, outputfilename):
        outputfile_path = Path(directory_path).joinpath(outputfilename)
        crash_file = Path(directory_path).joinpath("CRASH")
        if outputfile_path.exists():
            content = open(outputfile_path, "r").read()
            if ("Error" in content) or ("error" in content) or ("stopping" in content):
                logger.error(f"{outputfilename} output some wrong results. The program will exit!!!") 
                sys.exit(1)
            elif crash_file.exists():
                logger.error(f"{outputfilename} output some wrong results. The program will exit!!!") 
                sys.exit(1)
        else:
            logger.debug(f"{outputfilename} doesn't exist temporarily!")



