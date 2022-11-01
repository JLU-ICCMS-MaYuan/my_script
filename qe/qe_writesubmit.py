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
            npool=other_class.npool,
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
            npool=other_class.npool,
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
            npool=other_class.npool,

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
            npool=other_class.npool,
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
            npool=other_class.npool,
        )
        return self


    def write_submit_scripts(self, inpufilename, mode=None):

        if mode==None:
            mode=self.mode

        if mode == "relax-vc":
            jobname = self.s1_relax(self.work_underpressure, inpufilename)
            return jobname
        if mode == "scffit":
            jobname = self.s2_scffit(self.work_underpressure, inpufilename)
            return jobname
        if mode == "scf":
            jobname = self.s3_scf(self.work_underpressure, inpufilename)
            return jobname
        if mode =="nosplit":
            jobname = self.s4_PhNoSplit(self.work_underpressure, inpufilename)
            return jobname
        if mode =="split_dyn0":
            jobnames = []
            for i, inname in enumerate(inpufilename):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                jobname = self.s5_PhSplitDyn0(split_ph_dir, inname)
                jobnames.append(jobname)
                logger.info(f"finish submit job script in {i+1}")
            return jobnames
        if mode =="split_assignQ":
            jobnames = []
            for i, inname in enumerate(inpufilename):
                jobname = self.s5_PhSplitAssignQ(self.work_underpressure, inname)
                jobnames.append(jobname)
            return jobnames
        if mode =="q2r":
            jobname = self.s6_q2r(self.work_underpressure, inpufilename)
            return jobname
        if mode =="matdyn":
            jobname = self.s7_matdyn(self.work_underpressure, inpufilename)
            return jobname
        if mode =="matdyn_dos":
            jobname = self.s8_matdyn_dos(self.work_underpressure, inpufilename)
            return jobname
        if mode =="McAD":
            jobname = self.s9_lambda(self.work_underpressure, inpufilename)
            return jobname
        if mode =="eliashberg":
            jobname = self.s9_eliashberg(self.work_underpressure, inpufilename)
            return jobname
        if mode =="nscf":
            jobname = self.s10_nscf(self.work_underpressure, inpufilename)
            return jobname

    #  job scripts
    def s1_relax(self, _dirpath, inpufilename):
        _inpufilename = inpufilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s1_relax.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool {} <{}> {}  \n'.format(self.core, qebin_path,  self.npool, _inpufilename, _outputfilename))
            j.write('check symmetry ops is consistent or not after vc-relax                      \n')
            j.write('grep "Sym. Ops." relax.out                                                  \n')
            j.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")
        return jobname

    def s2_scffit(self, _dirpath, inpufilename):
        _inpufilename = inpufilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s2_scffit.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool {} <{}> {}  \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))                                                        
        return jobname
        
    def s3_scf(self, _dirpath, inpufilename):
        _inpufilename = inpufilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s3_scf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x -npool {} <{}> {} \n'.format(self.core, qebin_path, self.npool, _inpufilename,  _outputfilename)) 
        return jobname

    def s4_PhNoSplit(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s4_PhNoSplit.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/ph.x -npool {} <{}> {} \n'.format(self.core, qebin_path,  self.npool, _inpufilename, _outputfilename))
        return jobname

    def s5_PhSplitDyn0(self, _dirpath, inputfilename):
        _inputscffit_name, _inputscf_name, _inputsplitph_name = inputfilename
        _outputscffit_name  = _inputscffit_name.split(".")[0] + ".out"
        _outputscf_name     = _inputscf_name.split(".")[0] + ".out"
        _outputsplitph_name = _inputsplitph_name.split(".")[0] + ".out"
        jobname = "s5_PhSplitDyn0.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('echo "run scf.fit"                                                     \n')
            j.write('mpirun -np {} {}/pw.x -npool {} <{}> {} \n'.format(self.core, qebin_path,  self.npool, _inputscffit_name, _outputscffit_name))
            j.write('echo "run scf"                                                         \n')
            j.write('mpirun -np {} {}/pw.x -npool {} <{}> {} \n'.format(self.core, qebin_path,  self.npool, _inputscf_name,     _outputscf_name))
            j.write('echo "run split_ph"                                                    \n')
            j.write('mpirun -np {} {}/ph.x -npool {} <{}> {} \n'.format(self.core, qebin_path,  self.npool, _inputsplitph_name,  _outputsplitph_name))   
        return jobname

    def s5_PhSplitAssignQ(self, _dirpath, inputfilename):
        _inputsplitph_name = inputfilename
        _outputsplitph_name = _inputsplitph_name.split(".")[0] + ".out"
        jobname = "s5_"+_inputsplitph_name.split(".")[0]+".sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/ph.x -npool {} <{}> {}  \n'.format(self.core, qebin_path, self.npool, _inputsplitph_name, _outputsplitph_name))
        return jobname

    def s6_q2r(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s6_q2r.sh"
        _script_filepath = os.path.join(_dirpath,jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/q2r.x -npool {} <{}> {} \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))
            j.write('grep nqs q2r.out > nqs                   \n')  
        return jobname

    def s7_matdyn(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s7_matdyn.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/matdyn.x -npool {} <{}> {} \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))
        return jobname

    def s8_matdyn_dos(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s8_matdyn_dos.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/matdyn.x -npool {} <{}> {}  \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))
        return jobname

    def s9_lambda(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s9_lambda.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/lambda.x -npool {} <{}> {}  \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))
        return jobname

    def s9_eliashberg(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s9_eliashberg.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('killall -9 pw.x                   \n')
            j.write('\n\n                              \n')
            j.write('time {} > eliashberg.log 2>&1     \n'.format(eliashberg_x_path))
        return jobname

    def s10_nscf(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s10_nscf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('mpirun -np {} {}/pw.x {} <{}> {}  \n'.format(self.core, qebin_path, self.npool, _inpufilename, _outputfilename))
        return jobname

