import os

from qe_inputpara import qe_inputpara

class qe_submitjob:
    
    def __init__(self, qe_input_object, submit_job_system="slurm", run_mode=None):

        if isinstance(qe_input_object, qe_inputpara):
            self._qe_inputpara = qe_input_object
        self.submit_job_system = submit_job_system 
        self.run_mode          = run_mode

        if self.submit_job_system == "slurm":
            self.slurm_submitjob()
        elif self.submit_job_system == "pbs":
            self.pbs_submitjob()

    def pbs_job_system(self):
        if self.run_mode == "relax":
            self.pbsrelax(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scffit":
            self.pbsscfFit(self.submit_path)
        if self.run_mode == "scf":
            self.pbsscf(self.submit_path)
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
