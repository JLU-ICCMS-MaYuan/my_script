import re
import os
import sys
import time
import logging


from pathlib import Path

from qe.qe_inputpara import qe_inputpara
from qe.qebin import qebin_path, eliashberg_x_path

logger = logging.getLogger("qe_submitjob")

class qe_submitjob:
    
    def __init__(
        self,
        work_path: Path,
        submit_job_system: str,
        mode: str, 
        **kwargs
        ):

        self.work_path = work_path
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
            work_path=other_class.work_path,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
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
            work_path=other_class.work_path,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self 

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_path=other_class.work_path,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )

        return self 

    def submit_mode0(self, inputfilename, dotx_file):
        input_file = Path(self.work_path).joinpath(inputfilename)
        if not input_file.exists():
            raise FileExistsError(f" {inputfilename} doesn't exist")
        outputfilename = inputfilename.split(".")[0] + ".out"
        print("Note: --------------------")
        print("    !!!!!!!! Please Attention, You have been source your Intel Compiler !!!!!!!!")
        if dotx_file == "ph.x":
            cwd_path = os.getcwd()
            os.chdir(self.work_path)
            jobids = os.popen(f"nohup {qebin_path}/{dotx_file} <{inputfilename}> {outputfilename} 2>&1 & echo $!").read()
            print(f"{dotx_file} is running. pid or jobids = {jobids}")
            os.chdir(cwd_path)
            return jobids, outputfilename
        elif dotx_file == "lambda.x":
            cwd_path = os.getcwd()
            os.chdir(self.work_path)
            os.system(f"{qebin_path}/{dotx_file} <{inputfilename}> {outputfilename}")
            os.chdir(cwd_path)
        elif dotx_file == "eliashberg.x":
            cwd_path = os.getcwd()
            os.chdir(self.work_path)
            os.system(f"{eliashberg_x_path} > eliashberg.log")
            os.chdir(cwd_path)
        

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
        input_file = Path(self.work_path).joinpath(inputfilename)
        if not input_file.exists():
            raise FileExistsError(f" {inputfilename} doesn't exist")
        job_file = Path(self.work_path).joinpath(jobname)
        if not job_file.exists():
            raise FileExistsError(f" {jobname} doesn't exist")
        cwd = input_file.cwd()
        dst_dir = input_file.parent.absolute()
        os.chdir(dst_dir)
        if self.submit_job_system == "bash":
            print(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &")
            res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
            jobids = self.getpid()
        else:
            print(f"{self.submit_order} {jobname}")
            res = os.popen(f"{self.submit_order} {jobname}").read()
            jobids = re.findall(r"\d+", res)

        print(f"{jobname} is running. pid or jobids = {jobids}")
        os.chdir(cwd)
        # 检查任务是否成功提交，成功提交的话，应该会有进程号或者任务号返回。
        # 如果没有成功提交任务就跳出程序
        if not jobids:
            raise ValueError(f"The *.x of qe didn't run ! Because the jobids={jobids}. The program will exit! The order you use is {self.submit_order} {jobname}")
        return jobids

    def submit_mode2(self, inputfilename, jobname):
        """
        submit_mode1 can be used to submit:
           nosplit dyn0_flag=True     : only to get *dyn0 file. When *dyn0 appear, The ph.x will be kill  
           nosplit dyn0_flag=False    : run phono calculation, and don't interupte the ph.x
        """
        # 检查是否只是为了获得dyn0文件
        if self.dyn0_flag:
            # 检查是否只是为了获得dyn0文件，如果是：
            print(f"You set dyn0_flag = {self.dyn0_flag}. So When the dyn0 appears, the program will kill ph.x")
            # 在运行ph.x之前，检查dyn0文件是否已经存在
            if self.checksuffix(self.work_path, ".dyn0"):
                # 在运行ph.x之前，如果dyn0文件已经存在， 那么就直接退出程序
                logger.warning("Before running the ph.x, the program will check the *.dyn0 exists whether or not !")
                logger.warning("It seems that dyn0 is not create by you!! Please check it carefully!!! The program will exit!!!")
                sys.exit(0)
            else:
                # 在运行ph.x之前，如果dyn0文件不存在， 那么就执行ph.x的运行
                jobids, outputfilename = self.submit_mode0(inputfilename, "ph.x")
                # 如果执行完ph.x的运行后，检查返回的任务号不为空，说明ph.x的运行没有问题。
                while True:
                    # 然后检查dyn0文件是否存在，一旦产生就退出ph.x的运行。
                    time.sleep(3)
                    if self.checksuffix(self.work_path, ".dyn0"):
                        print("The *.dyn0 has been created just now !!! The program will run `killall -9 ph.x`")
                        os.system("killall -9 ph.x")
                        break
                    if self.checkerror(self.work_path, outputfilename):
                        print("The {outputfilename} has ERROR !!! The program will exit")
                        sys.exit(1)
        else:
            # 检查是否只是为了获得dyn0文件，如果否：
            print(f"You set dyn0_flag = {self.dyn0_flag}. The program will just simply run ph.x ")
            jobids = self.submit_mode1(inputfilename, jobname)

    def submit_mode3(self, inputfilename, jobnames):
        """split_dyn0"""
        jobids = []
        for i, jobname in enumerate(jobnames):
            cwd = os.getcwd()
            os.chdir(self.work_path.joinpath(str(i+1)))
            if self.submit_job_system == "bash":
                print(f"nohup {self.submit_order} {jobname} > phbash{i+1}.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                print(f"{self.submit_order} {jobname}")
                res = os.popen(f"{self.submit_order} {jobname}").read()
                jobids = re.findall(r"\d+", res)
            print(f"finish submit {jobname}, jobids = {''.join(jobids)}")
            os.chdir(cwd)

    def submit_mode4(self, inputfilename, jobnames):
        """split_assignQ"""
        jobids = []
        cwd = os.getcwd()
        os.chdir(self.work_path)
        for i, jobname in enumerate(jobnames):
            if self.submit_job_system == "bash":
                print(f"nohup {self.submit_order} {jobname} > phbash{i+1}.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                print(f"{self.submit_order} {jobname}")
                res = os.popen(f"{self.submit_order} {jobname}").read()
                jobids = re.findall(r"\d+", res)
                print(f"finish submit {jobname}, jobids = {''.join(jobids)}")
        os.chdir(cwd)



    @staticmethod
    def getpid():
        """get pid number"""
        jobids = []
        print("wait 6s, The program will tell you PID"); time.sleep(6); 
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
                print(f"{outputfilename} output some wrong results. The program will exit!!!") 
                sys.exit(1)
            elif crash_file.exists():
                print(f"{outputfilename} output some wrong results. The program will exit!!!") 
                sys.exit(1)
        else:
            print(f"{outputfilename} doesn't exist temporarily!")



