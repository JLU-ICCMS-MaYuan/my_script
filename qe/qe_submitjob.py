import re
import os
import logging
import time

from pathlib import Path
from unittest.util import sorted_list_difference
from qe_inputpara import qe_inputpara

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


    def submit(self, jobname):

        if self.mode == "relax-vc":
            dst_file = Path(self.work_underpressure).joinpath("relax.in")
            if not dst_file.exists():
                raise FileExistsError(" relax.in doesn't exist")

        cwd = dst_file.cwd()
        dst_dir = dst_file.parent.absolute()
        os.chdir(dst_dir)
        if self.submit_job_system == "bash":
            res   = os.popen(f"nohup {self.submit_order} {jobname} > bash.log 2>&1 &").read()
            jobid = os.popen(f"ps -aux | grep pw.x | awk '{print $2}'").read().split("\n"); print(jobid)
        else:
            res = os.popen(f"{self.submit_order} {jobname}").read()
            jobid = re.search(r"\d+", res)

        logger.info(f"pid or jobid = {jobid}")
        os.chdir(cwd)

        if self.mode == "scffit":
            dst_files = Path(self.work_underpressure).glob("scf.fit.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe scffit is running")
                    os.chdir(cwd)
        if self.mode == "scf":
            dst_files = Path(self.work_underpressure).glob("scf.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe scf is running")
                    os.chdir(cwd)
        if self.mode == "nscf":
            dst_files = Path(self.work_underpressure).glob("nscf.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe nscf is running")
                    os.chdir(cwd)
        if self.mode =="nosplit":
            dst_files = Path(self.work_underpressure).glob("ph_no_split.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    id_info = os.popen(f"{self.submit_order} {jobname}").read()
                    print(f"{id_info}")
                    if re.search(r"\d+", id_info) is not None:
                        id_num  = re.search(r"\d+", id_info).group()
                    os.chdir(cwd)
            while self.dyn0_flag:
                if os.path.exists(os.path.join(self.work_underpressure, self.system_name+".dyn0")):
                    os.system("sq")
                    os.system("scancel {}".format(id_num))
                    logger.info(f"The script detected the *.dyn0, so scancel the slurm job {id_num}")
                    self.dyn0_flag = False
                else:
                    time.sleep(5)
                    self.dyn0_flag = True
        if self.mode =="split_from_dyn0":
            for root, dirs, files in os.walk(self.work_underpressure):
                if "split_ph.in" in files and "scf.fit.in" in files and "scf.in" in files and "slurmph_split_from_dyn0.sh" in files:
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system(f"{self.submit_order} {jobname}")
                    os.chdir(cwd)
        if self.mode =="split_specify_q":
            split_ph_files = list(Path(self.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.qirreduced:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    cwd = os.getcwd()
                    os.chdir(self.work_underpressure)
                    os.system(f"{self.submit_order} {jobname}".format(split_ph_name))
                    logger.info(f"finish submit {split_ph_name}")
                    os.chdir(cwd) 
        if self.mode =="q2r":
            dst_files = Path(self.work_underpressure).glob("q2r.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system(f"{self.submit_order} {jobname}")
                    logger.info("qe q2r is running")
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
            dst_files = Path(self.work_underpressure).glob("matdyn.dos.in")
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
