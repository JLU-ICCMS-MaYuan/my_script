import re
import os
import logging
import time

from pathlib import Path

from qe.qe_inputpara import qe_inputpara

logger = logging.getLogger("qe_submitjob")

class qe_submitjob:
    
    def __init__(
        self,
        work_underpressure: Path,
        submit_job_system: str,
        mode: str, 
        **kwargs
        ):

        self.work_underpressure = work_underpressure
        self.submit_job_system  = submit_job_system
        self.mode = mode

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.submit_job_system == "slurm":
            self.submit_order = "sbatch"
        elif self.submit_job_system == "pbs":
            self.submit_order = "qsub"
        elif self.submit_job_system == "bash":
            self.submit_order = "bash"
        else:
            self.submit_order = ''


    @classmethod
    def init_from_relaxinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            dyn0_flag=other_class.dyn0_flag,
            system_name=other_class.system_name,
            qirreduced=other_class.qirreduced,
        )

        return self 

    @classmethod
    def init_from_dosinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self 

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self 


    def submit_mode1(self, inputfilename, jobname):
        """
        submit_mode1 can be used to submit:
            relax-vc, 
            scffit, 
            scf, 
            nscf,
            q2r,
            matdyn
        """
        input_file = Path(self.work_underpressure).joinpath(inputfilename)
        if not input_file.exists():
            raise FileExistsError(f" {inputfilename} doesn't exist")
        job_file = Path(self.work_underpressure).joinpath(jobname)
        if not job_file.exists():
            raise FileExistsError(f" {jobname} doesn't exist")
        cwd = input_file.cwd()
        dst_dir = input_file.parent.absolute()
        os.chdir(dst_dir)
        if self.submit_job_system == "bash":
            logger.info(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &")
            res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
            jobid = self.getpid()
        else:
            logger.info(f"{self.submit_order} {jobname}")
            res = os.popen(f"{self.submit_order} {jobname}").read()
            jobid = re.findall(r"\d+", res)
        logger.info(f"{jobname} is running. pid or jobid = {jobid}")
        os.chdir(cwd)

        return jobid

    def submit_mode2(self, inputfilename, jobname):
        """
        submit_mode1 can be used to submit:
           nosplit dyn0_flag=True     : only to get *dyn0 file. When *dyn0 appear, The ph.x will be kill  
           nosplit dyn0_flag=False    : run phono calculation, and don't interupte the ph.x
        """
        self.submit_mode1(inputfilename, jobname)
        logger.info(f"You set dyn0_flag = {self.dyn0_flag}.")
        if self.dyn0_flag:
            while True:
                time.sleep(8)
                if self.checksuffix(self.work_underpressure, ".dyn0"):
                    logger.info("The *.dyn0 has been created just now !!! The program will run `killall -9 ph.x`")
                    os.system("killall -9 ph.x")
                    break
                else:
                    logger.info("The *.dyn0 has existed ! It seems that dyn0 is not create by you!! Please check it carefully!!! The program will run `killall -9 ph.x`")
                    os.system("killall -9 ph.x")



    def submit_mode3(self, inputfilename, jobnames):
        """split_dyn0"""
        jobids = []
        for i, jobname in enumerate(jobnames):
            cwd = os.getcwd()
            os.chdir(self.work_underpressure.joinpath(str(i+1)))
            if self.submit_job_system == "bash":
                print(f"nohup {self.submit_order} {jobname} > phbash{i+1}.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                print(f"{self.submit_order} {jobname}")
                res = os.popen(f"{self.submit_order} {jobname}").read()
                jobids = re.findall(r"\d+", res)
            logger.info(f"finish submit {jobname}, jobid = {''.join(jobids)}")
            os.chdir(cwd)



    def submit_mode4(self, inputfilename, jobnames):
        """split_assignQ"""
        jobids = []
        cwd = os.getcwd()
        os.chdir(self.work_underpressure)
        for i, jobname in enumerate(jobnames):
            if self.submit_job_system == "bash":
                print(f"nohup {self.submit_order} {jobname} > phbash{i+1}.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                print(f"{self.submit_order} {jobname}")
                res = os.popen(f"{self.submit_order} {jobname}").read()
                jobids = re.findall(r"\d+", res)
                logger.info(f"finish submit {jobname}, jobid = {''.join(jobids)}")
        os.chdir(cwd)

        if self.mode =="matdyn":
            dst_files = Path(self.work_underpressure).glob("matdyn.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe matdyn is running")
                    os.chdir(cwd)
        if self.mode =="matdyn_dos":
            dst_files = Path(self.work_underpressure).glob("matdyn_dos.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe matdyn_dos is running")
                    os.chdir(cwd)
        if self.mode =="McAD":
            dst_files = Path(self.work_underpressure).glob("lambda.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe lambda is running")
                    os.chdir(cwd)
        if self.mode =="eliashberg":   
            dst_files = list(Path(self.work_underpressure).glob("ALPHA2F.OUT"))
            if len(dst_files) == 1:
                cwd = dst_files[0].cwd()
                dst_dir = dst_files[0].parent.absolute()
                os.chdir(dst_dir)
                os.system(f"{self.submit_order} {jobname}")
                logger.info("qe lambda is running")
                os.chdir(cwd)



    @staticmethod
    def getpid():
        """get pid number"""
        logger.info("wait 2s, The program will tell you PID"); time.sleep(2); 
        osawk = """ps -ef | grep -E "pw.x|ph.x|matdyn.x|lambda.x|q2r.x" |  grep -v grep | awk '{print $2}'""" # return a series of number, such as: 423423 324233 423424
        # ps -ef ps -ef用于查看全格式的全部进程，其中“ps”是在Linux中是查看进程的命令，“-e ”参数代表显示所有进程，“-f”参数代表全格式。
        # grep -E  ‘grep’ ‘-E’ 选项表示使用扩展的正则表达式。如果你使用 ‘grep’ 命令时带 ‘-E’，你只需要用途 ‘|’ 来分隔OR条件。 grep -E 'pattern1|pattern2' filename
        # grep -v grep 这里可以比较看出，多出了一个进程号，这是grep时所多出来的进程，通过grep -v grep去除包含grep文本的进程行 ，避免影响最终数据的正确性
        #  awk '{print $2}' 这样，就可以抓取PID号
        _jobid = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
        jobid = _jobid.strip("\n").split("\n")
        return jobid

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

