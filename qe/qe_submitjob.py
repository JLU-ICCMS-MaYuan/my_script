import re
import os
import logging
import time

from pathlib import Path
from qe_inputpara import qe_inputpara

logger = logging.getLogger("qe_workflow")

class qe_submitjob:
    
    def __init__(self, qe_input_object, submit_job_system="slurm", run_mode=None, **kwargs):

        if isinstance(qe_input_object, qe_inputpara):
            self._qe_inputpara = qe_input_object
        self.submit_job_system = submit_job_system 
        self.run_mode          = run_mode
        self.q_non_irreducible_amount = None
        if kwargs:
            for key, value in kwargs.items():
                if key == "dyn0_flag":
                    self.dyn0_flag = value
                if key == "q_non_irreducible_amount":
                    self.q_non_irreducible_amount = value
        
        if self.submit_job_system == "slurm":
            self.slurm_submitjob()
        elif self.submit_job_system == "pbs":
            self.pbs_submitjob()


    def slurm_submitjob(self):
        if self.run_mode == "relax":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("relax.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmrelax.sh")
                    logger.info("qe relax is running")
                    os.chdir(cwd)
        if self.run_mode == "scffit":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("scf.fit.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmscffit.sh")
                    logger.info("qe scffit is running")
                    os.chdir(cwd)
        if self.run_mode == "scf":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("scf.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmscf.sh")
                    logger.info("qe scf is running")
                    os.chdir(cwd)
        if self.run_mode =="ph_no_split":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("ph_no_split.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    id_info = os.popen("sbatch slurmph_no_split.sh").read()
                    print(f"{id_info}")
                    if re.search(r"\d+", id_info) is not None:
                        id_num  = re.search(r"\d+", id_info).group()
                    os.chdir(cwd)
            while self.dyn0_flag:
                if os.path.exists(os.path.join(self._qe_inputpara.work_underpressure, self._qe_inputpara.system_name+".dyn0")):
                    os.system("sq")
                    os.system("scancel {}".format(id_num))
                    logger.info(f"The script detected the *.dyn0, so scancel the slurm job {id_num}")
                    self.dyn0_flag = False
                else:
                    time.sleep(5)
                    self.dyn0_flag = True
        if self.run_mode =="ph_split_form_dyn0":
            for root, dirs, files in os.walk(self._qe_inputpara.work_underpressure):
                if "split_ph.in" in files and "scf.fit.in" in files and "scf.in" in files and "slurmph_split_form_dyn0.sh" in files:
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system("sbatch slurmph_split_form_dyn0.sh")
                    os.chdir(cwd)
        if self.run_mode =="ph_split_set_startlast_q":
            split_ph_files = list(Path(self._qe_inputpara.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.q_non_irreducible_amount:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    cwd = os.getcwd()
                    os.chdir(self._qe_inputpara.work_underpressure)
                    os.system("sbatch slurm_{}.sh".format(split_ph_name))
                    logger.info(f"finish submit {split_ph_name}")
                    os.chdir(cwd) 
        if self.run_mode =="q2r":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("q2r.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmq2r.sh")
                    logger.info("qe q2r is running")
                    os.chdir(cwd)
        if self.run_mode =="matdyn":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("matdyn.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmmatdyn.sh")
                    logger.info("qe matdyn is running")
                    os.chdir(cwd)
        if self.run_mode =="matdyn_dos":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("matdyn.dos.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmmatdyn_dos.sh")
                    logger.info("qe matdyn_dos is running")
                    os.chdir(cwd)
        if self.run_mode =="lambda":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("lambda.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("sbatch slurmlambda.sh")
                    logger.info("qe lambda is running")
                    os.chdir(cwd)

    def pbs_submitjob(self):
        if self.run_mode == "relax":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("relax.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsrelax.sh")
                    logger.info("qe relax is running")
                    os.chdir(cwd)
        if self.run_mode == "scffit":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("scf.fit.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsscffit.sh")
                    logger.info("qe scffit is running")
                    os.chdir(cwd)
        if self.run_mode == "scf":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("scf.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsscf.sh")
                    logger.info("qe scf is running")
                    os.chdir(cwd)
        if self.run_mode =="ph_no_split":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("ph_no_split.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    id_info = os.popen("qsub pbsph_no_split.sh").read()
                    print(f"{id_info}")
                    if re.search(r"\d+", id_info) is not None:
                        id_num  = re.search(r"\d+", id_info).group()
                    os.chdir(cwd)
            while self.dyn0_flag:
                if os.path.exists(os.path.join(self._qe_inputpara.work_underpressure, self._qe_inputpara.system_name+".dyn0")):
                    os.system("qstat")
                    os.system("qdel {}".format(id_num))
                    logger.info(f"The script detected the *.dyn0, so scancel the slurm job {id_num}")
                    self.dyn0_flag = False
                else:
                    time.sleep(5)
                    self.dyn0_flag = True
        if self.run_mode =="ph_split_form_dyn0":
            for root, dirs, files in os.walk(self._qe_inputpara.work_underpressure):
                if "split_ph.in" in files and "scf.fit.in" in files and "scf.in" in files and "pbsph_split_form_dyn0.sh" in files:
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system("qsub pbsph_split_form_dyn0.sh")
                    os.chdir(cwd)
        if self.run_mode =="ph_split_set_startlast_q":
            split_ph_files = list(Path(self._qe_inputpara.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.q_non_irreducible_amount:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    cwd = os.getcwd()
                    os.chdir(self._qe_inputpara.work_underpressure)
                    os.system("qsub pbs_{}.sh".format(split_ph_name))
                    logger.info(f"finish submit {split_ph_name}")
                    os.chdir(cwd) 
        if self.run_mode =="q2r":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("q2r.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsq2r.sh")
                    logger.info("qe q2r is running")
                    os.chdir(cwd)
        if self.run_mode =="matdyn":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("matdyn.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsmatdyn.sh")
                    logger.info("qe matdyn is running")
                    os.chdir(cwd)
        if self.run_mode =="matdyn_dos":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("matdyn.dos.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbsmatdyn_dos.sh")
                    logger.info("qe matdyn_dos is running")
                    os.chdir(cwd)
        if self.run_mode =="lambda":
            dst_files = Path(self._qe_inputpara.work_underpressure).glob("lambda.in")
            for dst_file in dst_files:
                if dst_file.exists():
                    cwd = dst_file.cwd()
                    dst_dir = dst_file.parent.absolute()
                    os.chdir(dst_dir)
                    os.system("qsub pbslambda.sh")
                    logger.info("qe lambda is running")
                    os.chdir(cwd)
