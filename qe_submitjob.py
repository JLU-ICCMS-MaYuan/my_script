import os
import logging

from pathlib import Path
from qe_inputpara import qe_inputpara

logger = logging.getLogger("qe_workflow")

class qe_submitjob:
    
    def __init__(self, qe_input_object, submit_job_system="slurm", run_mode=None):

        if isinstance(qe_input_object, qe_inputpara):
            self._qe_inputpara = qe_input_object
        self.submit_job_system = submit_job_system 
        self.run_mode          = run_mode
        self.mode2input_dict   = {
            "relax"  : "relax.in",
            "scffit" : "scf.fit.in",
            
        }

        if self.submit_job_system == "slurm":
            self.slurm_submitjob()
        elif self.submit_job_system == "pbs":
            self.pbs_submitjob()

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
            self.pbsph_no_split(self.submit_path)
        if self.run_mode =="ph_split_form_dyn0":
            self.pbsph_split_form_dyn0(self.submit_path)
        if self.run_mode =="ph_split_set_startlast_q":
            self.pbsph_split_set_startlast_q(self.submit_path, self.system_name)
        if self.run_mode =="q2r":
            self.pbsq2r(self.submit_path)
        if self.run_mode =="matdyn":
            self.pbsmatgen(self.submit_path)
        if self.run_mode =="matdyn_dos":
            self.pbsmatgen_dos(self.submit_path)
        if self.run_mode =="lambda":
            self.pbslambda(self.submit_path)
