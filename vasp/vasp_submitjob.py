import os
import re
import time
import shutil
import logging

from pathlib import Path

from vasp.vasp_inputpara import vasp_inputpara 

logger = logging.getLogger(__name__)

class vasp_submitjob:
    
    def __init__(
        self,
        vasp_inputpara:vasp_inputpara,
        ):

        self.vasp_inputpara = vasp_inputpara

        if self.vasp_inputpara.submit_job_system == "slurm":
            self.submit_order = "sbatch"
        elif self.vasp_inputpara.submit_job_system == "pbs":
            self.submit_order = "qsub"
        elif self.vasp_inputpara.submit_job_system == "lsf":
            self.submit_order = "bsub <"
        elif self.vasp_inputpara.submit_job_system == "bash":
            self.submit_order = "bash"
        else:
            self.submit_order = ''


    def submit_mode1(self, jobname, submit_path=None):
        """
        submit_mode1 can be used to submit:
            rvf
            rv3
            dfpt
            disp
        """
        if submit_path is None:
            submit_path = self.vasp_inputpara.work_path
            
        inputfilename = ["POSCAR", "POTCAR"]
        for input_name in inputfilename:
            input_file =submit_path.joinpath(input_name)
            if not input_file.exists():
                raise FileExistsError(f" {input_name} doesn't exist")
                
        job_file = submit_path.joinpath(jobname)
        if not job_file.exists():
            raise FileExistsError(f" {jobname} doesn't exist")
        cwd = input_file.cwd()
        dst_dir = input_file.parent.absolute()
        os.chdir(dst_dir)
        if self.vasp_inputpara.submit_job_system == "bash":
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
            raise ValueError(f"The vasp didn't run ! Because the jobids={jobids}. The program will exit! The order you use is {self.submit_order} {jobname}")
        return jobids

    def submit_mode2(self, jobname):

        patter = re.compile(r"POSCAR\-[0-9]{3}")
        poscar_files = os.listdir(self.vasp_inputpara.work_path)
        poscar_number_list = [patter.match(x).group() for x in poscar_files if patter.match(x)]
        for poscar_number in poscar_number_list:
            dst_number_dir = Path(self.vasp_inputpara.work_path).joinpath("disp-" + poscar_number.split("-")[-1])
            if not os.path.exists(dst_number_dir):
                os.makedirs(dst_number_dir)
            src_poscar = Path(self.vasp_inputpara.work_path).joinpath(poscar_number) ; dst_poscar = Path(dst_number_dir).joinpath("POSCAR"); shutil.copy(src_poscar, dst_poscar)
            src_potcar = Path(self.vasp_inputpara.work_path).joinpath("POTCAR")      ; dst_potcar = Path(dst_number_dir).joinpath("POTCAR"); shutil.copy(src_potcar, dst_potcar)
            src_incar  = Path(self.vasp_inputpara.work_path).joinpath("INCAR")       ; dst_incar  = Path(dst_number_dir).joinpath("INCAR" ); shutil.copy(src_incar, dst_incar )
            src_kpoints= Path(self.vasp_inputpara.work_path).joinpath("KPOINTS")     ; dst_kpoints= Path(dst_number_dir).joinpath("KPOINTS"); shutil.copy(src_kpoints, dst_kpoints) if src_kpoints.exists() else None
            src_submit = Path(self.vasp_inputpara.work_path).joinpath(jobname)       ; dst_submit = Path(dst_number_dir).joinpath(jobname);  shutil.copy(src_submit, dst_submit)
            cwd = src_poscar.cwd()
            dst_dir = dst_number_dir.absolute()
            os.chdir(dst_dir)
            if self.vasp_inputpara.submit_job_system == "bash":
                print(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &")
                res = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
                jobids = self.getpid()
            else:
                print(f"{self.submit_order} {jobname}")
                res = os.popen(f"{self.submit_order} {jobname}").read()
                jobids = re.findall(r"\d+", res)

            print(f"{jobname} is running. pid or jobids = {jobids}")
            os.chdir(cwd)


    @staticmethod
    def getpid():
        """get pid number"""
        jobids = []
        print("wait 3s, The program will tell you PID"); time.sleep(3); 
        osawk = """ps -ef | grep -E "vasp_std" |  grep -v grep | awk '{print $2}'""" # return a series of number, such as: 423423 324233 423424
        # ps -ef ps -ef用于查看全格式的全部进程，其中“ps”是在Linux中是查看进程的命令，“-e ”参数代表显示所有进程，“-f”参数代表全格式。
        # grep -E  ‘grep’ ‘-E’ 选项表示使用扩展的正则表达式。如果你使用 ‘grep’ 命令时带 ‘-E’，你只需要用途 ‘|’ 来分隔OR条件。 grep -E 'pattern1|pattern2' filename
        # grep -v grep 这里可以比较看出，多出了一个进程号，这是grep时所多出来的进程，通过grep -v grep去除包含grep文本的进程行 ，避免影响最终数据的正确性
        #  awk '{print $2}' 这样，就可以抓取PID号
        _jobids = os.popen(osawk).read()  # return a string; such as '423423\n324233\n423424\n'
        jobids = _jobids.strip("\n").split("\n")
        return jobids