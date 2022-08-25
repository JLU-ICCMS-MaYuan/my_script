from asyncio.log import logger
import os
import re
import shutil
from vasp_inputpara import vasp_inputpara 


class vasp_submitjob:
    
    def __init__(
        self,
        work_underpressure: str, 
        submit_job_system: str, 
        mode: str,
        ):

        self.work_underpressure = work_underpressure
        self.submit_job_system = submit_job_system
        self.mode = mode

        if self.submit_job_system == "slurm":
            self.slurm_submitjob()
        elif self.submit_job_system == "pbs":
            self.pbs_submitjob()

    @classmethod
    def init_from_relaxinput(cls, other_class: vasp_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )
        
        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: vasp_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
        )
        
        return self

    def slurm_submitjob(self):
        if self.mode == "rvf":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("sbatch slurmFopt.sh")
            os.chdir(cwd)
        elif self.mode == "rv3":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("sbatch slurm3opt.sh")
            os.chdir(cwd)
        elif self.mode == "disp":
            patter = re.compile(r"POSCAR\-[0-9]{3}")
            poscar_files = os.listdir(self.work_underpressure)
            poscar_number_list = [patter.match(x).group() for x in poscar_files if patter.match(x)]
            for poscar_number in poscar_number_list:
                dst_number_dir = os.path.join(self.work_underpressure, "disp-" + poscar_number.split("-")[-1])
                if not os.path.exists(dst_number_dir):
                    os.makedirs(dst_number_dir)
                src_poscar = os.path.join(self.work_underpressure, poscar_number) ; dst_poscar = os.path.join(dst_number_dir, "POSCAR");     shutil.copy(src_poscar, dst_poscar)
                src_potcar = os.path.join(self.work_underpressure, "POTCAR")      ; dst_potcar = os.path.join(dst_number_dir, "POTCAR");     shutil.copy(src_potcar, dst_potcar)
                src_incar  = os.path.join(self.work_underpressure, "INCAR_disp")  ; dst_incar  = os.path.join(dst_number_dir, "INCAR" );     shutil.copy(src_incar, dst_incar )
                src_kpoints= os.path.join(self.work_underpressure, "KPOINTS")     ; dst_kpoints= os.path.join(dst_number_dir,"KPOINTS");     shutil.copy(src_kpoints, dst_kpoints)
                src_slurm  = os.path.join(self.work_underpressure, "slurmdisp.sh"); dst_slurm  = os.path.join(dst_number_dir,"slurmdisp.sh");shutil.copy(src_slurm, dst_slurm)
                cwd = os.getcwd()
                os.chdir(dst_number_dir)
                os.system("sbatch slurmdisp.sh")
                os.chdir(cwd) 
        elif self.mode == "dfpt":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("sbatch slurmdfpt.sh")
            os.chdir(cwd) 

    def pbs_submitjob(self):
        if self.mode == "rvf":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("qsub pbsFopt.sh")
            logger.info(" vasp optfine is running.")
            os.chdir(cwd)
        elif self.mode == "rv3":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("qsub pbs3opt.sh")
            os.chdir(cwd)
            logger.info(" vasp opt3 is running.")
        elif self.mode == "disp":
            patter = re.compile(r"POSCAR\-[0-9]{3}")
            poscar_files = os.listdir(self.work_underpressure)
            poscar_number_list = [patter.match(x).group() for x in poscar_files if patter.match(x)]
            for poscar_number in poscar_number_list:
                dst_number_dir = os.path.join(self.work_underpressure, "disp-" + poscar_number.split("-")[-1])
                if not os.path.exists(dst_number_dir):
                    os.makedirs(dst_number_dir)
                src_poscar = os.path.join(self.work_underpressure, poscar_number) ; dst_poscar = os.path.join(dst_number_dir, "POSCAR");     shutil.copy(src_poscar, dst_poscar)
                src_potcar = os.path.join(self.work_underpressure, "POTCAR")      ; dst_potcar = os.path.join(dst_number_dir, "POTCAR");     shutil.copy(src_potcar, dst_potcar)
                src_incar  = os.path.join(self.work_underpressure, "INCAR_disp")  ; dst_incar  = os.path.join(dst_number_dir, "INCAR" );     shutil.copy(src_incar, dst_incar )
                src_kpoints= os.path.join(self.work_underpressure, "KPOINTS")     ; dst_kpoints= os.path.join(dst_number_dir,"KPOINTS");     shutil.copy(src_kpoints, dst_kpoints)
                src_pbs    = os.path.join(self.work_underpressure, "pbsdisp.sh")  ; dst_pbs  = os.path.join(dst_number_dir,"pbsdisp.sh");    shutil.copy(src_pbs, dst_pbs)
                cwd = os.getcwd()
                os.chdir(dst_number_dir)
                os.system("qsub pbsdisp.sh")
                logger.info(" vasp phono-disp is running.")
                os.chdir(cwd) 
        elif self.mode == "dfpt":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system("qsub pbsdfpt.sh")
            logger.info(" vasp phono-dfpt is running.")
            os.chdir(cwd) 