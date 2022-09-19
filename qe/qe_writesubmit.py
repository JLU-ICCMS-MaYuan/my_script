import os
import re
import logging
from pathlib import Path
from qe_inputpara import qe_inputpara
from qebin import qebin_path, eliashberg_x_path

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
            self.slurm_job_system()
        elif self.submit_job_system == "pbs":
            self.pbs_job_system()


    @classmethod
    def init_from_relaxinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
        )
        return self

    @classmethod
    def init_from_scfinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
        )
        return self
    
    @classmethod
    def init_from_phonoinput(cls, other_class: qe_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            qirreduced=other_class.qirreduced,
            qirreduced_coords=other_class.qirreduced_coords,

        )
        return self

    @classmethod
    def init_from_scinput(cls, other_class: qe_inputpara):

        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
        )
        return self



    def slurm_job_system(self):
        if self.mode == "relax-vc":
            self.slurmrelax(self.work_underpressure)
        if self.mode == "scffit":
            self.slurmscffit(self.work_underpressure)
        if self.mode == "scf":
            self.slurmscf(self.work_underpressure)
        if self.mode =="nosplit":
            self.slurmph_no_split(self.work_underpressure)
        if self.mode =="split_from_dyn0":
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                self.slurmph_split_from_dyn0(split_ph_dir)
                logger.info(f"finish submit job script in {i+1}")
        if self.mode =="split_specify_q":
            split_ph_files = list(Path(self.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.qirreduced:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    self.slurmph_split_set_startlast_q(self.work_underpressure, split_ph_name)
        if self.mode =="q2r":
            self.slurmq2r(self.work_underpressure)
        if self.mode =="matdyn":
            self.slurmmatdyn(self.work_underpressure)
        if self.mode =="matdyn_dos":
            self.slurmmatdyn_dos(self.work_underpressure)
        if self.mode =="McAD":
            self.slurmlambda(self.work_underpressure)
        if self.mode =="eliashberg":
            self.slurmeliashberg(self.work_underpressure)

    def pbs_job_system(self):
        if self.mode == "relax-vc":
            self.pbsrelax(self.work_underpressure)
        if self.mode == "scffit":
            self.pbsscffit(self.work_underpressure)
        if self.mode == "scf":
            self.pbsscf(self.work_underpressure)
        if self.mode =="nosplit":
            self.pbsph_no_split(self.work_underpressure)
        if self.mode =="split_from_dyn0":
            for i, q3 in enumerate(self.qirreduced_coords):
                split_ph_dir = os.path.join(self.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                self.pbsph_split_from_dyn0(split_ph_dir)
                logger.info(f"finish submit job script in {i+1}")
        if self.mode =="split_specify_q":
            split_ph_files = list(Path(self.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.qirreduced:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    self.pbsph_split_set_startlast_q(self.work_underpressure, split_ph_name)
        if self.mode =="q2r":
            self.pbsq2r(self.work_underpressure)
        if self.mode =="matdyn":
            self.pbsmatdyn(self.work_underpressure)
        if self.mode =="matdyn_dos":
            self.pbsmatdyn_dos(self.work_underpressure)
        if self.mode =="McAD":
            self.pbslambda(self.work_underpressure)
        if self.mode =="eliashberg":
            self.pbseliashberg(self.work_underpressure)



    # slurm job scripts
    def slurmrelax(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmrelax.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=relax                                                \n')                         
            slurm.write('#SBATCH  --output=log.relax.out                                          \n')                       
            slurm.write('#SBATCH  --error=log.relax.err                                           \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <relax.in> relax.out                        \n'.format(qebin_path))
            slurm.write('check symmetry ops is consistent or not after vc-relax                   \n')
            slurm.write('grep "Sym. Ops." relax.out                                               \n')
            slurm.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")

    def slurmscffit(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmscffit.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=scf.fit                                              \n')                         
            slurm.write('#SBATCH  --output=log.scf.fit.out                                        \n')                       
            slurm.write('#SBATCH  --error=log.scf.fit.err                                         \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <scf.fit.in> scf.fit.out                    \n'.format(qebin_path))                                                                         

    def slurmscf(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmscf.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=scf                                                  \n')                         
            slurm.write('#SBATCH  --output=log.scf.out                                            \n')                       
            slurm.write('#SBATCH  --error=log.scf.err                                             \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <scf.in> scf.out                            \n'.format(qebin_path))   

    def slurmnscf(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmnscf.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=nscf                                                 \n')                         
            slurm.write('#SBATCH  --output=log.nscf.out                                           \n')                       
            slurm.write('#SBATCH  --error=log.nscf.err                                            \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <nscf.in> nscf.out                          \n'.format(qebin_path))   

    def slurmph_no_split(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmph_no_split.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=ph_no_split                                          \n')                         
            slurm.write('#SBATCH  --output=log.ph_no_split.out                                    \n')                       
            slurm.write('#SBATCH  --error=log.ph_no_split.err                                     \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}ph.x -npool 4 <ph_no_split.in> ph_no_split.out            \n'.format(qebin_path))

    def slurmph_split_from_dyn0(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmph_split_from_dyn0.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=ph_split                                             \n')                         
            slurm.write('#SBATCH  --output=log.ph_split.out                                       \n')                       
            slurm.write('#SBATCH  --error=log.ph_split.err                                        \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('echo "run scf.fit"                                                       \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <scf.fit.in> scf.fit.out                    \n'.format(qebin_path))
            slurm.write('echo "run scf"                                                           \n')
            slurm.write('mpirun -n 48 {}pw.x -npool 4 <scf.in> scf.out                            \n'.format(qebin_path))
            slurm.write('echo "run split_ph"                                                      \n')
            slurm.write('mpirun -n 48 {}ph.x -npool 4 <split_ph.in> split_ph.out                  \n'.format(qebin_path))   

    def slurmph_split_set_startlast_q(self, slurm_dirpath, split_ph_name):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurm_"+split_ph_name+".sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                 \n')     
            slurm.write('#SBATCH  --job-name={}                                                    \n'.format(split_ph_name))                         
            slurm.write('#SBATCH  --output=log.{}.out                                              \n'.format(split_ph_name))                       
            slurm.write('#SBATCH  --error=log.{}.err                                               \n'.format(split_ph_name))                      
            slurm.write('#SBATCH  --partition={}                                                   \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                        \n')             
            slurm.write('#SBATCH  --ntasks=48                                                      \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                             \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                                \n')                     
            slurm.write('\n\n                                                                      \n')
            slurm.write('source /work/env/intel2018                                                \n')
            slurm.write('ulimit -s unlimited                                                       \n')
            slurm.write('\n\n                                                                      \n')
            # TODO 
            slurm.write('mpirun -n 48 {}ph.x -npool 4 <{}.in> {}.out                               \n'.format(qebin_path ,split_ph_name, split_ph_name))

    def slurmq2r(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmq2r.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=q2r                                                  \n')                         
            slurm.write('#SBATCH  --output=log.q2r.out                                            \n')                       
            slurm.write('#SBATCH  --error=log.q2r.err                                             \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}q2r.x -npool 4 <q2r.in> q2r.out                           \n'.format(qebin_path))
            slurm.write('grep nqs q2r.out > nqs                                                   \n')  

    def slurmmatdyn(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmmatdyn.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=matgen                                               \n')                         
            slurm.write('#SBATCH  --output=log.matgen.out                                         \n')                       
            slurm.write('#SBATCH  --error=log.matgen.err                                          \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}matdyn.x -npool 4 <matdyn.in> matdyn.out                  \n'.format(qebin_path))  

    def slurmmatdyn_dos(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmmatdyn_dos.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=matgen_dos                                           \n')                         
            slurm.write('#SBATCH  --output=log.matgen_dos.out                                     \n')                       
            slurm.write('#SBATCH  --error=log.matgen_dos.err                                      \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out          \n'.format(qebin_path))  

    def slurmlambda(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmlambda.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=lambda                                               \n')                         
            slurm.write('#SBATCH  --output=log.lambda.out                                         \n')                       
            slurm.write('#SBATCH  --error=log.lambda.err                                          \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 {}lambda.x <lambda.in> lambda.out                           \n'.format(qebin_path))  

    def slurmeliashberg(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "pbseliashberg.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=eliashberg                                           \n')                         
            slurm.write('#SBATCH  --output=log.eliashberg.out                                     \n')                       
            slurm.write('#SBATCH  --error=log.eliashberg.err                                      \n')                      
            slurm.write('#SBATCH  --partition={}                                                  \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write(' source /work/home/mayuan/intel/oneapi/setvars.sh --force                \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('cd $PBS_O_WORKDIR                                                        \n')
            slurm.write('killall -9 pw.x                                                          \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('time {} > eliashberg.log 2>&1                                            \n'.format(eliashberg_x_path))  

        


    # pbs job scripts
    def pbsrelax(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsrelax.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    relax                                                         \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -np 28 {}pw.x -npool 4 <relax.in> relax.out                       \n'.format(qebin_path))
            pbs.write('check symmetry ops is consistent or not after vc-relax                   \n')
            pbs.write('grep "Sym. Ops." relax.out                                               \n')
            pbs.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")

    def pbsscffit(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsscffit.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    scffit                                                        \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -np 28 {}pw.x -npool 4 <scf.fit.in> scf.fit.out                   \n'.format(qebin_path))

    def pbsscf(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsscf.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    scf                                                           \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}pw.x -npool 4 <scf.in> scf.out                            \n'.format(qebin_path))   

    def pbsnscf(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsnscf.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    nscf                                                          \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}pw.x -npool 4 <nscf.in> nscf.out                          \n'.format(qebin_path))   

    def pbsph_no_split(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsph_no_split.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    ph_no_split                                                   \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}ph.x -npool 4 <ph_no_split.in> ph_no_split.out            \n'.format(qebin_path))

    def pbsph_split_from_dyn0(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsph_split_from_dyn0.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    split_ph                                                      \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('echo "run scf.fit"                                                       \n')
            pbs.write('mpirun -n 28 {}pw.x -npool 4 <scf.fit.in> scf.fit.out                    \n'.format(qebin_path))
            pbs.write('echo "run scf"                                                           \n')
            pbs.write('mpirun -n 28 {}pw.x -npool 4 <scf.in> scf.out                            \n'.format(qebin_path))
            pbs.write('echo "run split_ph"                                                      \n')
            pbs.write('mpirun -n 28 {}ph.x -npool 4 <split_ph.in> split_ph.out                  \n'.format(qebin_path))   

    def pbsph_split_set_startlast_q(self, pbs_dirpath, split_ph_name):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbs_"+split_ph_name+".sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    {}                                                            \n'.format(split_ph_name))                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            # TODO 
            pbs.write('mpirun -n 28 {}ph.x -npool 4 <{}.in> {}.out                              \n'.format(qebin_path, split_ph_name, split_ph_name))

    def pbsq2r(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsq2r.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    q2r                                                           \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}q2r.x -npool 4 <q2r.in> q2r.out                           \n'.format(qebin_path))
            pbs.write('grep nqs q2r.out > nqs                                                   \n')  

    def pbsmatdyn(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsmatdyn.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    matdyn                                                        \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}matdyn.x -npool 4 <matdyn.in> matdyn.out                  \n'.format(qebin_path))  

    def pbsmatdyn_dos(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsmatdyn_dos.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    matdyn_dos                                                    \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out          \n'.format(qebin_path))  

    def pbslambda(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbslambda.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    lambda                                                        \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 28 {}lambda.x <lambda.in> lambda.out                           \n'.format(qebin_path))  

    def pbseliashberg(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbseliashberg.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    eliashberg                                                    \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('cd $PBS_O_WORKDIR                                                        \n')
            pbs.write('killall -9 pw.x                                                          \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('time {} > eliashberg.log 2>&1                                            \n'.format(eliashberg_x_path))  
        

