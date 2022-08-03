import os
import re
import logging
from pathlib import Path
from qe_inputpara import qe_inputpara

logger = logging.getLogger("qe_writesubmit")

class qe_writesubmit:

    def __init__(self, qe_input_object, submit_job_system="slurm", run_mode=None, **kwargs):

        if isinstance(qe_input_object, qe_inputpara):
            self._qe_inputpara  = qe_input_object
        self.submit_job_system        = submit_job_system 
        self.run_mode                 = run_mode
        self.q_non_irreducible_amount = None

        if kwargs:
            for key, value in kwargs.items():
                if key=="q_non_irreducible_amount":
                    self.q_non_irreducible_amount = value

        if self.submit_job_system == "slurm":
            self.slurm_job_system()
        elif self.submit_job_system == "pbs":
            self.pbs_job_system()

    def slurm_job_system(self):
        if self.run_mode == "relax":
            self.slurmrelax(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scffit":
            self.slurmscffit(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scf":
            self.slurmscf(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_no_split":
            self.slurmph_no_split(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_split_form_dyn0":
            dyn0_names = list(Path(self._qe_inputpara.work_underpressure).glob("*.dyn0"))
            if len(dyn0_names)==1:
                dyn0_path = str(dyn0_names[0].absolute())
            else:
                raise FileExistsError("Exist many *.dyn0 files or No *.dyn0")
            _, _, q_coordinate_list, _ = self._qe_inputpara.get_q_from_dyn0(dyn0_path)
            for i, q3 in enumerate(q_coordinate_list):
                split_ph_dir = os.path.join(self._qe_inputpara.work_underpressure, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                self.slurmph_split_form_dyn0(split_ph_dir)
                logger.info(f"finish submit job script in {i+1}")
        if self.run_mode =="ph_split_set_startlast_q":
            split_ph_files = list(Path(self._qe_inputpara.work_underpressure).glob("split_ph*.in"))
            if len(split_ph_files)==self.q_non_irreducible_amount:
                for split_ph_file in split_ph_files:
                    split_ph_name = re.split(r"[\/.]" ,str(split_ph_file))[-2]
                    self.slurmph_split_set_startlast_q(self._qe_inputpara.work_underpressure, split_ph_name)
        if self.run_mode =="q2r":
            self.slurmq2r(self._qe_inputpara.work_underpressure)
        if self.run_mode =="matdyn":
            self.slurmmatgen(self._qe_inputpara.work_underpressure)
        if self.run_mode =="matdyn_dos":
            self.slurmmatgen_dos(self._qe_inputpara.work_underpressure)
        if self.run_mode =="lambda":
            self.slurmlambda(self._qe_inputpara.work_underpressure)

    def pbs_job_system(self):
        if self.run_mode == "relax":
            self.pbsrelax(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scffit":
            self.pbsscffit(self._qe_inputpara.work_underpressure)
        if self.run_mode == "scf":
            self.pbsscf(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_no_split":
            self.pbsph_no_split(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_split_form_dyn0":
            self.pbsph_split_form_dyn0(self._qe_inputpara.work_underpressure)
        if self.run_mode =="ph_split_set_startlast_q":
            self.pbsph_split_set_startlast_q(self._qe_inputpara.work_underpressure, self._qe_inputpara.system_name)
        if self.run_mode =="q2r":
            self.pbsq2r(self._qe_inputpara.work_underpressure)
        if self.run_mode =="matdyn":
            self.pbsmatgen(self._qe_inputpara.work_underpressure)
        if self.run_mode =="matdyn_dos":
            self.pbsmatgen_dos(self._qe_inputpara.work_underpressure)
        if self.run_mode =="lambda":
            self.pbslambda(self._qe_inputpara.work_underpressure)


    # slurm job scripts
    def slurmrelax(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmrelax.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=relax                                             \n')                         
            slurm.write('#SBATCH  --output=log.relax.out                                      \n')                       
            slurm.write('#SBATCH  --error=log.relax.err                                       \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <relax.in> relax.out \n')
            slurm.write('check symmetry ops is consistent or not after vc-relax                   \n')
            slurm.write('grep "Sym. Ops." relax.out                                               \n')
            slurm.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")

    def slurmscffit(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmscffit.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=scf.fit                                             \n')                         
            slurm.write('#SBATCH  --output=log.scf.fit.out                                      \n')                       
            slurm.write('#SBATCH  --error=log.scf.fit.err                                       \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out \n')                                                                         

    def slurmscf(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmscf.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=scf                                                  \n')                         
            slurm.write('#SBATCH  --output=log.scf.out                                           \n')                       
            slurm.write('#SBATCH  --error=log.scf.err                                            \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out\n')   

    def slurmph_no_split(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmph_no_split.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                                 \n')     
            slurm.write('#SBATCH  --job-name=ph_no_split                                                           \n')                         
            slurm.write('#SBATCH  --output=log.ph_no_split.out                                                     \n')                       
            slurm.write('#SBATCH  --error=log.ph_no_split.err                                                      \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                                                \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                                        \n')             
            slurm.write('#SBATCH  --ntasks=48                                                                      \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                                             \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                                                \n')                     
            slurm.write('\n\n                                                                                      \n')
            slurm.write('source /work/env/intel2018                                                                \n')
            slurm.write('ulimit -s unlimited                                                                       \n')
            slurm.write('\n\n                                                                                      \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <ph_no_split.in> ph_no_split.out \n')

    def slurmph_split_form_dyn0(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmph_split_form_dyn0.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                          \n')     
            slurm.write('#SBATCH  --job-name=ph_split                                                       \n')                         
            slurm.write('#SBATCH  --output=log.ph_split.out                                                 \n')                       
            slurm.write('#SBATCH  --error=log.ph_split.err                                                  \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                                         \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                                 \n')             
            slurm.write('#SBATCH  --ntasks=48                                                               \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                                      \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                                         \n')                     
            slurm.write('\n\n                                                                               \n')
            slurm.write('source /work/env/intel2018                                                         \n')
            slurm.write('ulimit -s unlimited                                                                \n')
            slurm.write('\n\n                                                                               \n')
            slurm.write('echo "run scf.fit"                                                                 \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out  \n')
            slurm.write('echo "run scf"                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out          \n')
            slurm.write('echo "run split_ph"                                                                \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <split_ph.in> split_ph.out\n')   

    def slurmph_split_set_startlast_q(self, slurm_dirpath, split_ph_name):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurm_"+split_ph_name+".sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                          \n')     
            slurm.write('#SBATCH  --job-name={}                                                             \n'.format(split_ph_name))                         
            slurm.write('#SBATCH  --output=log.{}.out                                                       \n'.format(split_ph_name))                       
            slurm.write('#SBATCH  --error=log.{}.err                                                        \n'.format(split_ph_name))                      
            slurm.write('#SBATCH  --partition=xieyu                                                         \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                                 \n')             
            slurm.write('#SBATCH  --ntasks=48                                                               \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                                      \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                                         \n')                     
            slurm.write('\n\n                                                                               \n')
            slurm.write('source /work/env/intel2018                                                         \n')
            slurm.write('ulimit -s unlimited                                                                \n')
            slurm.write('\n\n                                                                               \n')
            # TODO 
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <{}.in> {}.out            \n'.format(split_ph_name, split_ph_name))

    def slurmq2r(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmq2r.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                  \n')     
            slurm.write('#SBATCH  --job-name=q2r                                                    \n')                         
            slurm.write('#SBATCH  --output=log.q2r.out                                         \n')                       
            slurm.write('#SBATCH  --error=log.q2r.err                                          \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                                 \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                         \n')             
            slurm.write('#SBATCH  --ntasks=48                                                       \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                              \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                                 \n')                     
            slurm.write('\n\n                                                                       \n')
            slurm.write('source /work/env/intel2018                                                 \n')
            slurm.write('ulimit -s unlimited                                                        \n')
            slurm.write('\n\n                                                                       \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/q2r.x -npool 4 <q2r.in> q2r.out \n')
            slurm.write('grep nqs q2r.out > nqs                                                     \n')  

    def slurmmatgen(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmmatgen.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=matgen                                               \n')                         
            slurm.write('#SBATCH  --output=log.matgen.out                                    \n')                       
            slurm.write('#SBATCH  --error=log.matgen.err                                     \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/matdyn.x -npool 4 <matdyn.in> matdyn.out \n')  

    def slurmmatgen_dos(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmmatgen_dos.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=matgen_dos                                           \n')                         
            slurm.write('#SBATCH  --output=log.matgen_dos.out                                    \n')                       
            slurm.write('#SBATCH  --error=log.matgen_dos.err                                     \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out \n')  

    def slurmlambda(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmlambda.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                                                \n')     
            slurm.write('#SBATCH  --job-name=lambda                                           \n')                         
            slurm.write('#SBATCH  --output=log.lambda.out                                    \n')                       
            slurm.write('#SBATCH  --error=log.lambda.err                                     \n')                      
            slurm.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                                                       \n')             
            slurm.write('#SBATCH  --ntasks=48                                                     \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            slurm.write('\n\n                                                                     \n')
            slurm.write('source /work/env/intel2018                                               \n')
            slurm.write('ulimit -s unlimited                                                      \n')
            slurm.write('\n\n                                                                     \n')
            slurm.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/lambda.x <lambda.in> lambda.out \n')  


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
            pbs.write('mpirun -np 28 /public/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <relax.in> relax.out \n')
            pbs.write('check symmetry ops is consistent or not after vc-relax                   \n')
            pbs.write('grep "Sym. Ops." relax.out                                               \n')
            pbs.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")

    def pbsscffit(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsscffit.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#PBS -N    scffit                                                         \n')                         
            pbs.write('#PBS -q    liuhy                                                         \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                                \n')             
            pbs.write('#PBS -j    oe                                                            \n')                        
            pbs.write('#PBS -V                                                                  \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /public/home/mayuan/intel/oneapi/setvars.sh                       \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -np 28 /public/home/mayuan/software/qe-7.1/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out \n')

    def pbsscf(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsscf.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=scf                                                  \n')                         
            pbs.write('#SBATCH  --output=log.scf.out                                           \n')                       
            pbs.write('#SBATCH  --error=log.scf.err                                            \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out\n')   

    def pbsph_no_split(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsph_no_split.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                                 \n')     
            pbs.write('#SBATCH  --job-name=ph_no_split                                                           \n')                         
            pbs.write('#SBATCH  --output=log.ph_no_split.out                                                     \n')                       
            pbs.write('#SBATCH  --error=log.ph_no_split.err                                                      \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                                                \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                                        \n')             
            pbs.write('#SBATCH  --ntasks=48                                                                      \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                                             \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                                                \n')                     
            pbs.write('\n\n                                                                                      \n')
            pbs.write('source /work/env/intel2018                                                                \n')
            pbs.write('ulimit -s unlimited                                                                       \n')
            pbs.write('\n\n                                                                                      \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <ph_no_split.in> ph_no_split.out \n')

    def pbsph_split_form_dyn0(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsph_split_form_dyn0.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                          \n')     
            pbs.write('#SBATCH  --job-name=ph_split                                                       \n')                         
            pbs.write('#SBATCH  --output=log.ph_split.out                                                 \n')                       
            pbs.write('#SBATCH  --error=log.ph_split.err                                                  \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                                         \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                                 \n')             
            pbs.write('#SBATCH  --ntasks=48                                                               \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                                      \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                                         \n')                     
            pbs.write('\n\n                                                                               \n')
            pbs.write('source /work/env/intel2018                                                         \n')
            pbs.write('ulimit -s unlimited                                                                \n')
            pbs.write('\n\n                                                                               \n')
            pbs.write('echo "run scf.fit"                                                                 \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out  \n')
            pbs.write('echo "run scf"                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.in> scf.out          \n')
            pbs.write('echo "run split_ph"                                                                \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <split_ph.in> split_ph.out\n')   

    def pbsph_split_set_startlast_q(self, pbs_dirpath, split_ph_name):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbs_"+split_ph_name+".sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                          \n')     
            pbs.write('#SBATCH  --job-name={}                                                             \n'.format(split_ph_name))                         
            pbs.write('#SBATCH  --output=log.{}.out                                                       \n'.format(split_ph_name))                       
            pbs.write('#SBATCH  --error=log.{}.err                                                        \n'.format(split_ph_name))                      
            pbs.write('#SBATCH  --partition=xieyu                                                         \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                                 \n')             
            pbs.write('#SBATCH  --ntasks=48                                                               \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                                      \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                                         \n')                     
            pbs.write('\n\n                                                                               \n')
            pbs.write('source /work/env/intel2018                                                         \n')
            pbs.write('ulimit -s unlimited                                                                \n')
            pbs.write('\n\n                                                                               \n')
            # TODO 
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/ph.x -npool 4 <{}.in> {}.out            \n'.format(split_ph_name, split_ph_name))

    def pbsq2r(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsq2r.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                  \n')     
            pbs.write('#SBATCH  --job-name=q2r                                                    \n')                         
            pbs.write('#SBATCH  --output=log.q2r.out                                         \n')                       
            pbs.write('#SBATCH  --error=log.q2r.err                                          \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                                 \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                         \n')             
            pbs.write('#SBATCH  --ntasks=48                                                       \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                              \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                                 \n')                     
            pbs.write('\n\n                                                                       \n')
            pbs.write('source /work/env/intel2018                                                 \n')
            pbs.write('ulimit -s unlimited                                                        \n')
            pbs.write('\n\n                                                                       \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/q2r.x -npool 4 <q2r.in> q2r.out \n')
            pbs.write('grep nqs q2r.out > nqs                                                     \n')  

    def pbsmatgen(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsmatgen.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=matgen                                               \n')                         
            pbs.write('#SBATCH  --output=log.matgen.out                                    \n')                       
            pbs.write('#SBATCH  --error=log.matgen.err                                     \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/matdyn.x -npool 4 <matdyn.in> matdyn.out \n')  

    def pbsmatgen_dos(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsmatgen_dos.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=matgen_dos                                           \n')                         
            pbs.write('#SBATCH  --output=log.matgen_dos.out                                    \n')                       
            pbs.write('#SBATCH  --error=log.matgen_dos.err                                     \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/matdyn.x -npool 4 <matdyn.dos.in> matdyn.dos.out \n')  

    def pbslambda(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbslambda.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=lambda                                           \n')                         
            pbs.write('#SBATCH  --output=log.lambda.out                                    \n')                       
            pbs.write('#SBATCH  --error=log.lambda.err                                     \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/lambda.x <lambda.in> lambda.out \n')  
        

