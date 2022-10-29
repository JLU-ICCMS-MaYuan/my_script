import os
import re
import logging
from pathlib import Path
from qe_inputpara import qe_inputpara
from qebin import qebin_path, qe_source_libs, eliashberg_x_path, bashtitle, slurmtitle, pbstitle

logger = logging.getLogger("qe_writesubmit")

class qe_writesubmit:

    def __init__(
        self,
        work_underpressure: Path,
        submit_job_system: str,
        mode: str,
        **kwargs,
        ):

        self.work_underpressure = work_underpressure
        self.submit_job_system  = submit_job_system 
        self.mode               = mode

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.submit_job_system == "slurm":
            self.jobtitle =  slurmtitle, 
        elif self.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.submit_job_system == "bash":
            self.jobtitle = bashtitle
        else:
            self.jobtitle = ''

    @classmethod
    def init_from_relaxinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        return self
    
    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,

            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,
        )
        return self

    @classmethod
    def init_from_dosinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        return self

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):

        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        return self


    def job_system(self):

        if self.mode == "relax-vc":
            jobname = self.s1_relax(self.work_underpressure)
        if self.mode == "scffit":
            jobname = self.s2_scffit(self.work_underpressure)
        if self.mode == "scf":
            jobname = self.s3_scf(self.work_underpressure)
        if self.mode =="nosplit":
            jobname = self.s4_PhNoSplit(self.work_underpressure)
        if self.mode =="split_from_dyn0":
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                jobname = self.s5_PhSplitFromDyn0(split_ph_dir)
                logger.info(f"finish submit job script in {i+1}")
        if self.mode =="split_specify_q":
            split_ph_files = list(Path(self.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.qirreduced:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    jobname = self.s5_PhSplitSetStartLastQ(self.work_underpressure, split_ph_name)
        if self.mode =="q2r":
            jobname = self.s6_q2r(self.work_underpressure)
        if self.mode =="matdyn":
            jobname = self.s7_matdyn(self.work_underpressure)
        if self.mode =="matdyn_dos":
            jobname = self.s8_matdyn_dos(self.work_underpressure)
        if self.mode =="McAD":
            jobname = self.s9_lambda(self.work_underpressure)
        if self.mode =="eliashberg":
            jobname = self.s9_eliashberg(self.work_underpressure)
        if self.mode =="nscf":
            jobname = self.s10_nscf(self.work_underpressure)

        return jobname

    #  job scripts
    def s1_relax(self, _dirpath):
        jobname = "s1_relax.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool 4 <relax.in> relax.out                          \n'.format(self.core, qebin_path))
            j.write('check symmetry ops is consistent or not after vc-relax                      \n')
            j.write('grep "Sym. Ops." relax.out                                                  \n')
            j.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")
        return jobname

    def s2_scffit(self, _dirpath):
        jobname = "s2_scffit.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool 4 <scf.fit.in> scf.fit.out                    \n'.format(self.core, qebin_path))                                                                         
        return jobname
        
    def s3_scf(self, _dirpath):
        jobname = "s3_scf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool 4 <scf.in> scf.out                            \n'.format(self.core, qebin_path))   
        return jobname

    def s4_PhNoSplit(self, _dirpath):
        jobname = "s4_PhNoSplit.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/ph.x -npool 4 <ph_no_split.in> ph_no_split.out            \n'.format(self.core, qebin_path))
        return jobname

    def s5_PhSplitFromDyn0(self, _dirpath):
        jobname = "s5_PhSplitFromDyn0.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('echo "run scf.fit"                                                       \n')
            j.write('mpirun -np {} {}/pw.x -npool 4 <scf.fit.in> scf.fit.out                    \n'.format(self.core, qebin_path))
            j.write('echo "run scf"                                                           \n')
            j.write('mpirun -np {} {}/pw.x -npool 4 <scf.in> scf.out                            \n'.format(self.core, qebin_path))
            j.write('echo "run split_ph"                                                      \n')
            j.write('mpirun -np {} {}/ph.x -npool 4 <split_ph.in> split_ph.out                  \n'.format(self.core, qebin_path))   
        return jobname

    def s5_PhSplitSetStartLastQ(self, _dirpath, split_ph_name):
        jobname = split_ph_name+".sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/ph.x -npool 4 <{}.in> {}.out                               \n'.format(self.core, qebin_path ,split_ph_name, split_ph_name))
        return jobname

    def s6_q2r(self, _dirpath):
        jobname = "s6_q2r.sh"
        _script_filepath = os.path.join(_dirpath,jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/q2r.x -npool 4 <q2r.in> q2r.out                           \n'.format(self.core, qebin_path))
            j.write('grep nqs q2r.out > nqs                                                   \n')  
        return jobname

    def s7_matdyn(self, _dirpath):
        jobname = "s7_matdyn.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/matdyn.x -npool 4 <matdyn.in> matdyn.out                  \n'.format(self.core, qebin_path))  
        return jobname

    def s8_matdyn_dos(self, _dirpath):
        jobname = "s8_matdyn_dos.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out          \n'.format(self.core, qebin_path))  

    def s9_lambda(self, _dirpath):
        jobname = "s9_lambda.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/lambda.x <lambda.in> lambda.out                           \n'.format(self.core, qebin_path))  
        return jobname

    def s9_eliashberg(self, _dirpath):
        jobname = "s9_eliashberg.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('killall -9 pw.x                                                          \n')
            j.write('\n\n                                                                     \n')
            j.write('time {} > eliashberg.log 2>&1                                            \n'.format(self.core, eliashberg_x_path))  
        return jobname

    def s10_nscf(self, _dirpath):
        jobname = "s10_nscf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.write(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool 4 <nscf.in> nscf.out                          \n'.format(self.core, qebin_path))
        return jobname

