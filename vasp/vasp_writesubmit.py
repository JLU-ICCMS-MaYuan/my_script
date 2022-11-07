import os

from vasp.vasp_inputpara import vasp_inputpara
from vasp.vaspbin import vaspbin_path, bashtitle, slurmtitle, pbstitle


class vasp_writesubmit:
    
    def __init__(
        self, 
        work_underpressure: str, 
        submit_job_system: str, 
        mode: str,
        queue: str,
        core: int,
        ):
        self.work_underpressure = work_underpressure
        self.submit_job_system = submit_job_system
        self.mode = mode
        self.queue = queue
        self.core = core

        if self.submit_job_system == "slurm":
            self.jobtitle =  slurmtitle
        elif self.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.submit_job_system == "bash":
            self.jobtitle = bashtitle
        else:
            self.jobtitle = ''

    @classmethod
    def init_from_relaxinput(cls, other_class: vasp_inputpara):
        
        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        
        return self
    
    @classmethod
    def init_from_phonoinput(cls, other_class: vasp_inputpara):

        self = cls(
            work_underpressure=other_class.work_underpressure,
            submit_job_system=other_class.submit_job_system,
            mode=other_class.mode,
            queue=other_class.queue,
            core=other_class.core,
        )
        
        return self

    def write_submit_scripts(self, mode=None):

        if mode==None:
            mode=self.mode

        if mode == 'rvf':
            self.Fopt(self.work_underpressure)
        elif mode == 'rv3':
            self.threeopt(self.work_underpressure)
        elif mode == 'disp':
            self.disp(self.work_underpressure)
        elif mode == 'dfpt':
            self.dfpt(self.work_underpressure)


    # submit job scripts
    def Fopt(self, submit_dirpath):
        jobname = "Fopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('cp INCAR_fine INCAR                           \n')                    
            submit.write('num=0                                         \n')                                                                               
            submit.write('while true;do                                 \n')              
            submit.write('        let num+=1                            \n')                   
            submit.write('        echo "run fine vasp opt-$num"         \n')                                      
            submit.write('        killall -9 vasp_std                                \n')                                                                         
            submit.write('        sleep 3                                            \n')                                                                         
            submit.write('        timeout 14400s mpirun -n {} {} > vasp.log 2>&1    \n'.format(self.core, vaspbin_path))                                                                         
            submit.write('        cp -f CONTCAR CONTCAR-fine &&  cp -f CONTCAR POSCAR\n')                                                            
            submit.write("        rows=`sed -n '/F\=/p' OSZICAR | wc -l`             \n")                                               
            submit.write('        echo "rows-$rows"                                  \n')                           
            submit.write('        echo $num >> count_opt_times                       \n')
            submit.write('        if [ "$rows" -eq "1" ];then                        \n')                                    
            submit.write('                break                                      \n')                      
            submit.write('        fi                                                 \n')           
            submit.write('        if [ "$num" -gt "50" ]; then                       \n') 
            submit.write('                break                                      \n')       
            submit.write('        fi                                                 \n')        
            submit.write('done                                                       \n')  
        return jobname


    def threeopt(self, submit_dirpath):
        jobname = "threeopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("for i in {1..3}; do                  \n")
            submit.write("cp INCAR_$i INCAR                    \n")
            submit.write("killall -9 vasp_std                  \n")
            submit.write("sleep 3                              \n")
            submit.write("timeout 14400s mpirun -np {} {} > vasp.log_$i 2>&1\n".format(self.core, vaspbin_path))
            submit.write("cp CONTCAR POSCAR                    \n")
            submit.write("done                                 \n")
        return jobname
         

    def dfpt(self, submit_dirpath):
        jobname = "dfpt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('echo "run fine DFPT"                \n')
            submit.write('cp -f INCAR_dfpt INCAR              \n')
            submit.write('cp -f SPOSCAR POSCAR                \n')
            submit.write('mpirun -np {} {} > vasp.log 2>&1    \n'.format(self.core, vaspbin_path))                        
        return jobname


    def disp(self, submit_dirpath):
        jobname = "disp.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        print(submit_script_filepath)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('echo "run Displacement" && pwd      \n')
            submit.write('mpirun -np {} {} > vasp.log 2>&1  \n'.format(self.core, vaspbin_path))   
        return jobname
 