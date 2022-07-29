import os

    
class qe_submitjob:

    def __init__(self):
        pass

    @classmethod
    def pbsrelax(cls, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsrelax.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=relax                                             \n')                         
            pbs.write('#SBATCH  --output=log.relax.out                                      \n')                       
            pbs.write('#SBATCH  --error=log.relax.err                                       \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <relax.in> relax.out \n')
            pbs.write('check symmetry ops is consistent or not after vc-relax                   \n')
            pbs.write('grep "Sym. Ops." relax.out                                               \n')
            pbs.write("awk '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out \n")

    @classmethod
    def pbsscfFit(cls, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsscfFit.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                                \n')     
            pbs.write('#SBATCH  --job-name=scf.fit                                             \n')                         
            pbs.write('#SBATCH  --output=log.scf.fit.out                                      \n')                       
            pbs.write('#SBATCH  --error=log.scf.fit.err                                       \n')                      
            pbs.write('#SBATCH  --partition=xieyu                                               \n')    # lhy lbt is both ok                
            pbs.write('#SBATCH  --nodes=1                                                       \n')             
            pbs.write('#SBATCH  --ntasks=48                                                     \n')               
            pbs.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
            pbs.write('#SBATCH  --cpus-per-task=1                                               \n')                     
            pbs.write('\n\n                                                                     \n')
            pbs.write('source /work/env/intel2018                                               \n')
            pbs.write('ulimit -s unlimited                                                      \n')
            pbs.write('\n\n                                                                     \n')
            pbs.write('mpirun -n 48 /work/software/q-e-qe-6.8/bin/pw.x -npool 4 <scf.fit.in> scf.fit.out \n')                                                                         

    @classmethod
    def pbsscf(cls, pbs_dirpath):
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

    @classmethod
    def pbsph_no_split(cls, pbs_dirpath):
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

    @classmethod
    def pbsph_split_form_dyn0(cls, pbs_dirpath):
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

    @classmethod
    def pbsph_split_set_startlast_q(cls, pbs_dirpath, split_ph_name):
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

    @classmethod
    def pbsq2r(cls, pbs_dirpath):
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

    @classmethod
    def pbsmatgen(cls, pbs_dirpath):
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

    @classmethod
    def pbsmatgen_dos(cls, pbs_dirpath):
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

    @classmethod
    def pbslambda(cls, pbs_dirpath):
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
        