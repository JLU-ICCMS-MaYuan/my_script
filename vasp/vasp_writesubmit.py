import os

from vasp.vasp_inputpara import vasp_inputpara
from vasp.vaspbin import vaspbin_path, intel_compiler


class vasp_writesubmit:
    
    def __init__(
        self, 
        work_underpressure: str, 
        submit_job_system: str, 
        mode: str,
        queue: str,
        ):
        self.work_underpressure = work_underpressure
        self.submit_job_system = submit_job_system
        self.mode = mode
        self.queue = queue

        if self.submit_job_system == 'slurm':
            self.slurm_job_system()
        elif self.submit_job_system == 'pbs':
            self.pbs_job_system()

    @classmethod
    def init_from_relaxinput(cls, other_class: vasp_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
        )
        
        return self
    
    @classmethod
    def init_from_phonoinput(cls, other_class: vasp_inputpara):

        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
        )
        
        return self

    def slurm_job_system(self):

        if self.mode == 'rvf':
            self.slurmFopt(self.work_underpressure)
        elif self.mode == 'rv3':
            self.slurm3opt(self.work_underpressure)
        elif self.mode == 'disp':
            self.slurmdisp(self.work_underpressure)
        elif self.mode == 'dfpt':
            self.slurmdfpt(self.work_underpressure)

    def pbs_job_system(self):
        
        if self.mode == 'rvf':
            self.pbsFopt(self.work_underpressure)
        elif self.mode == 'rv3':
            self.pbs3opt(self.work_underpressure)
        elif self.mode == 'disp':
            self.pbsdisp(self.work_underpressure)
        elif self.mode == 'dfpt':
            self.pbsdfpt(self.work_underpressure)

    # slurm job scripts
    def slurmFopt(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmFopt.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write('#!/bin/sh                                     \n')     
            slurm.write('#SBATCH  --job-name=opt_fine                  \n')                         
            slurm.write('#SBATCH  --output=opt_fine.out.%j             \n')                       
            slurm.write('#SBATCH  --error=opt_fine.err.%j              \n')                      
            slurm.write('#SBATCH  --partition={}                       \n'.format(self.queue))    # lhy lbt is both ok                
            slurm.write('#SBATCH  --nodes=1                            \n')             
            slurm.write('#SBATCH  --ntasks=48                          \n')               
            slurm.write('#SBATCH  --ntasks-per-node=48                 \n')                        
            slurm.write('#SBATCH  --cpus-per-task=1                    \n')                     
            slurm.write('\n\n                                          \n')
            slurm.write('source {}                                     \n'.format(intel_compiler))
            slurm.write('ulimit -s unlimited                           \n')
            slurm.write('export I_MPI_ADJUST_REDUCE=3                  \n')
            slurm.write('export MPIR_CVAR_COLL_ALIAS_CHECK=0           \n')
            slurm.write('\n\n                                          \n')
            slurm.write('cp INCAR_fine INCAR                           \n')                    
            slurm.write('num=0                                         \n')                                                                               
            slurm.write('while true;do                                 \n')              
            slurm.write('        let num+=1                            \n')                   
            slurm.write('        echo "run fine vasp opt-$num"         \n')                                      
            slurm.write('        killall -9 vasp_std                                \n')                                                                         
            slurm.write('        sleep 3                                            \n')                                                                         
            slurm.write('        timeout 14400s mpirun -n 48 {} > vasp.log 2>&1    \n'.format(vaspbin_path))                                                                         
            slurm.write('        cp -f CONTCAR CONTCAR-fine &&  cp -f CONTCAR POSCAR\n')                                                            
            slurm.write("        rows=`sed -n '/F\=/p' OSZICAR | wc -l`             \n")                                               
            slurm.write('        echo "rows-$rows"                                  \n')                           
            slurm.write('        echo $num >> count_opt_times                       \n')
            slurm.write('        if [ "$rows" -eq "1" ];then                        \n')                                    
            slurm.write('                break                                      \n')                      
            slurm.write('        fi                                                 \n')           
            slurm.write('        if [ "$num" -gt "50" ]; then                       \n') 
            slurm.write('                break                                      \n')       
            slurm.write('        fi                                                 \n')        
            slurm.write('done                                                       \n')  

    def slurm3opt(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurm3opt.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write("#!/bin/sh                            \n")     
            slurm.write("#SBATCH  --job-name=opt3steps        \n")                         
            slurm.write("#SBATCH  --output=opt3steps.out.%j   \n")                       
            slurm.write("#SBATCH  --error=opt3steps.err.%j    \n")                      
            slurm.write("#SBATCH  --partition={}              \n".format(self.queue))                   
            slurm.write("#SBATCH  --nodes=1                   \n")             
            slurm.write("#SBATCH  --ntasks=48                 \n")               
            slurm.write("#SBATCH  --ntasks-per-node=48        \n")                        
            slurm.write("#SBATCH  --cpus-per-task=1           \n")                     
            slurm.write("\n                                   \n")
            slurm.write('source {}                            \n'.format(intel_compiler))

            slurm.write("ulimit -s unlimited                  \n")        
            slurm.write("export I_MPI_ADJUST_REDUCE=3         \n")        
            slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0  \n")

            slurm.write("                                     \n")
            slurm.write("for i in {1..3}; do                  \n")
            slurm.write("cp INCAR_$i INCAR                    \n")
            slurm.write("killall -9 vasp_std                  \n")
            slurm.write("sleep 3                              \n")
            slurm.write("timeout 14400s mpirun -np 48 {} > vasp.log_$i 2>&1\n".format(vaspbin_path))
            slurm.write("cp CONTCAR POSCAR                    \n")
            slurm.write("done                                 \n")
          
         

    def slurmdfpt(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmdfpt.sh")
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write("#!/bin/sh                           \n")     
            slurm.write("#SBATCH  --job-name=dfpt            \n")                         
            slurm.write("#SBATCH  --output=dfpt.out.%j       \n")                       
            slurm.write("#SBATCH  --error=dfpt.err.%j        \n")                      
            slurm.write("#SBATCH  --partition={}             \n".format(self.queue))    # lhy lbt is both ok                
            slurm.write("#SBATCH  --nodes=1                  \n")             
            slurm.write("#SBATCH  --ntasks=48                \n")               
            slurm.write("#SBATCH  --ntasks-per-node=48       \n")                        
            slurm.write("#SBATCH  --cpus-per-task=1          \n")                     
            slurm.write("                                    \n")
            slurm.write('source {}                           \n'.format(intel_compiler))
            slurm.write("ulimit -s unlimited                 \n")
            slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
            slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
            slurm.write("                                    \n")
            slurm.write('echo "run fine DFPT"                \n')
            slurm.write('cp -f INCAR_dfpt INCAR              \n')
            slurm.write('cp -f SPOSCAR POSCAR                \n')
            slurm.write('mpirun -np 48 {} > vasp.log 2>&1    \n'.format(vaspbin_path))                        

    def slurmdisp(self, slurm_dirpath):
        slurm_script_filepath = os.path.join(slurm_dirpath, "slurmdisp.sh")
        print(slurm_script_filepath)
        with open(slurm_script_filepath, "w") as slurm:
            slurm.write("#!/bin/sh                           \n")     
            slurm.write("#SBATCH  --job-name=disp            \n")                         
            slurm.write("#SBATCH  --output=disp.out.%j       \n")                       
            slurm.write("#SBATCH  --error=disp.err.%j        \n")                      
            slurm.write("#SBATCH  --partition={}             \n".format(self.queue))    # lhy lbt is both ok                
            slurm.write("#SBATCH  --nodes=1                  \n")             
            slurm.write("#SBATCH  --ntasks=48                \n")               
            slurm.write("#SBATCH  --ntasks-per-node=48       \n")                        
            slurm.write("#SBATCH  --cpus-per-task=1          \n")                     
            slurm.write("                                  \n\n")
            slurm.write('source {}                           \n'.format(intel_compiler))
            slurm.write("ulimit -s unlimited                 \n")
            slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
            slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
            slurm.write("                                  \n\n")
            slurm.write('echo "run Displacement" && pwd      \n')
            slurm.write('mpirun -np 48 {} > vasp.log 2>&1  \n'.format(vaspbin_path))   


    # pbs job scripts
    def pbsFopt(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsFopt.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                           \n')     
            pbs.write('#PBS -N    Fopt                                     \n')                         
            pbs.write('#PBS -q    liuhy                                    \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                           \n')             
            pbs.write('#PBS -j    oe                                       \n')                        
            pbs.write('#PBS -V                                             \n')                     
            pbs.write('\n\n                                                \n')
            pbs.write('source {}                                           \n'.format(intel_compiler))
            
            pbs.write('ulimit -s unlimited                                 \n')
            pbs.write("export I_MPI_ADJUST_REDUCE=3        \n")
            pbs.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")

            pbs.write('cd $PBS_O_WORKDIR                                          \n')
            pbs.write('killall -9 pw.x                                            \n')
            pbs.write('\n                                                         \n')
            pbs.write('cp INCAR_fine INCAR                                        \n')                    
            pbs.write('num=0                                                      \n')                                                                               
            pbs.write('while true;do                                              \n')              
            pbs.write('        let num+=1                                         \n')                   
            pbs.write('        echo "run fine vasp opt-$num"                      \n')                                      
            pbs.write("        killall -9 vasp_std                                \n")
            pbs.write("        sleep 3                                            \n")
            pbs.write('        timeout 14400s mpirun -n 28 {}  > vasp.log 2>&1    \n'.format(vaspbin_path))                        
            pbs.write('        cp -f CONTCAR CONTCAR-fine &&  cp -f CONTCAR POSCAR\n')                                                            
            pbs.write("        rows=`sed -n '/F\=/p' OSZICAR | wc -l`             \n")                                               
            pbs.write('        echo "rows-$rows"                                  \n')                           
            pbs.write('        echo $num >> count_opt_times                       \n')
            pbs.write('        if [ "$rows" -eq "1" ];then                        \n')                                    
            pbs.write('                break                                      \n')                      
            pbs.write('        fi                                                 \n')           
            pbs.write('        if [ "$num" -gt "50" ]; then                       \n') 
            pbs.write('                break                                      \n')       
            pbs.write('        fi                                                 \n')        
            pbs.write('done                                                       \n')  

    def pbs3opt(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbs3opt.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                  \n')     
            pbs.write('#PBS -N    3opt                                            \n')                         
            pbs.write('#PBS -q    liuhy                                           \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                  \n')             
            pbs.write('#PBS -j    oe                                              \n')                        
            pbs.write('#PBS -V                                                    \n')                     
            pbs.write('\n                                                         \n')
            pbs.write('source {}                                                  \n'.format(intel_compiler))
            
            pbs.write("ulimit -s unlimited                                        \n")
            pbs.write("export I_MPI_ADJUST_REDUCE=3                               \n")
            pbs.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0                        \n")

            pbs.write('cd $PBS_O_WORKDIR                                          \n')
            pbs.write('killall -9 pw.x                                            \n')
            pbs.write('\n                                                         \n')

            pbs.write("for i in {1..3}; do                                        \n")
            pbs.write("cp INCAR_$i INCAR                                          \n")
            pbs.write("killall -9 vasp_std                                        \n")
            pbs.write("sleep 3                                                    \n")
            pbs.write("timeout 14400s mpirun -np  28 {} > vasp.log_$i 2>&1      \n".format(vaspbin_path))
            pbs.write("cp CONTCAR POSCAR                                          \n")
            pbs.write("done                                                       \n")

    def pbsdfpt(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsdfpt.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                  \n')     
            pbs.write('#PBS -N    dfpt                                            \n')                         
            pbs.write('#PBS -q    liuhy                                           \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                  \n')             
            pbs.write('#PBS -j    oe                                              \n')                        
            pbs.write('#PBS -V                                                    \n')                     
            pbs.write('\n\n                                                       \n')
            pbs.write('source {}                                                  \n'.format(intel_compiler))
            
            pbs.write("ulimit -s unlimited                                        \n")
            pbs.write("export I_MPI_ADJUST_REDUCE=3                               \n")
            pbs.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0                        \n")

            pbs.write('cd $PBS_O_WORKDIR                                          \n')
            pbs.write('killall -9 pw.x                                            \n')
            pbs.write('\n                                                         \n')
            pbs.write('echo "run fine DFPT"                                       \n')
            pbs.write('cp -f INCAR_dfpt INCAR                                     \n')
            pbs.write('cp -f SPOSCAR POSCAR                                       \n')
            pbs.write('mpirun -n 28 {}  > vasp.log 2>&1                           \n'.format(vaspbin_path))                        

    def pbsdisp(self, pbs_dirpath):
        pbs_script_filepath = os.path.join(pbs_dirpath, "pbsdisp.sh")
        with open(pbs_script_filepath, "w") as pbs:
            pbs.write('#!/bin/sh                                                  \n')     
            pbs.write('#PBS -N    disp                                            \n')                         
            pbs.write('#PBS -q    liuhy                                           \n')    # lhy lbt is both ok                
            pbs.write('#PBS -l    nodes=1:ppn=28                                  \n')             
            pbs.write('#PBS -j    oe                                              \n')                        
            pbs.write('#PBS -V                                                    \n')                     
            pbs.write('\n                                                         \n')
            pbs.write('source {}                                                  \n'.format(intel_compiler))
            
            pbs.write("ulimit -s unlimited                                        \n")
            pbs.write("export I_MPI_ADJUST_REDUCE=3                               \n")
            pbs.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0                        \n")

            pbs.write('cd $PBS_O_WORKDIR                                          \n')
            pbs.write('killall -9 pw.x                                            \n')
            pbs.write('                                                           \n')
            pbs.write('echo "run Displacement" && pwd                             \n')
            pbs.write('mpirun -n 28 {}  > vasp.log 2>&1                           \n'.format(vaspbin_path))                        





