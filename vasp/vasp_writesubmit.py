import os
import re

from vasp.vasp_inputpara import vasp_inputpara
from vasp.vaspbin import vaspstd_path, vaspgam_path, bashtitle, slurmtitle, pbstitle


class vasp_writesubmit:
    
    def __init__(
        self, 
        vasp_inputpara
        ):

        self.vasp_inputpara = vasp_inputpara

        if self.vasp_inputpara.submit_job_system == "slurm":
            self.jobtitle = self.update_slurmPartition(slurmtitle)
        elif self.vasp_inputpara.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.vasp_inputpara.submit_job_system == "bash":
            self.jobtitle = bashtitle
        else:
            self.jobtitle = ''

    def update_slurmPartition(self, title:str, new_partition:str):
        if self.check_partition_exists(new_partition):
            # 使用正则表达式替换 --partition 后面的内容
            updated_title = re.sub(r'--partition=\S+', f'--partition={new_partition}', title)
            print(f"Partition exist! {new_partition}")
            return updated_title
        else:
            print(f"{new_partition} doesn't exist! Keep partition name in ~/.my_scripts.py")
            return title
        
    def check_partition_exists(self, new_partition: str) -> bool:
        """检查队列是否存在"""
        try:
            # 使用 sinfo 命令检查队列是否存在
            partitions = os.popen('sinfo -h --format=%P').read().splitlines()
            
            # 判断队列是否在返回的队列列表中
            if new_partition in partitions:
                return True
            else:
                return False
        except Exception as e:
            print(f"Error throws up when run `sinfo -h --format=%P`: {e}")
            return False

    def write_submit_scripts(self, mode=None, submitjob_path=None):

        if mode == None:
            mode = self.vasp_inputpara.mode
        if submitjob_path == None:
            submitjob_path = self.vasp_inputpara.work_path

        if mode == 'rvf':
            jobname = self.fopt(submitjob_path)
            return jobname
        elif mode == "rv1":
            jobname = self.oneopt(submitjob_path)
            return jobname
        elif mode == 'rv3':
            jobname = self.threeopt(submitjob_path)
            return jobname
        elif mode == 'rv4':
            jobname = self.fouropt(submitjob_path)
            return jobname
        elif mode == 'disp':
            jobname = self.disp(submitjob_path)
            return jobname
        elif mode == 'dfpt':
            jobname = self.dfpt(submitjob_path)
            return jobname
        elif mode in ['only-scf', "only-eband", 'only-eledos', 'only-cohp']:
            jobname = self.only_onemode(submitjob_path)
            return jobname
        elif mode == "scf-eband-eledos":
            jobname = self.scf_eband_eledos(submitjob_path)
            return jobname
        elif mode == "scf-eband":
            jobname = self.scf_eband(submitjob_path)
            return jobname
        elif mode == "scf-eledos":
            jobname = self.scf_eledos(submitjob_path)
            return jobname
        elif mode == "eband-eledos":
            jobname = self.eband_eledos(submitjob_path)
            return jobname
        elif mode == 'nvt':
            jobname = self.md_nvt(submitjob_path)
            return jobname
        elif mode == 'npt':
            jobname = self.md_npt(submitjob_path)
            return jobname
        elif mode == 'nve':
            jobname = self.md_nve(submitjob_path)
            return jobname
        
    # submit job scripts
    def fopt(self, submit_dirpath):
        jobname = "fopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('num=0                                         \n')                                                                               
            submit.write('while true;do                                 \n')              
            submit.write('        let num+=1                            \n')                   
            submit.write('        echo "run fine vasp opt-$num"         \n')                                      
            submit.write('        killall -9 vasp_std                                \n')                                                                         
            submit.write('        sleep 3                                            \n')                                                                         
            submit.write('        timeout 14400s {} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                                                                         
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

    def oneopt(self, submit_dirpath):
        jobname = "oneopt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                                                                         
            submit.write('cp CONTCAR POSCAR\n')                                                                
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
            submit.write("{} {} > vasp.log_$i 2>&1             \n".format(self.vasp_inputpara.execmd, vaspstd_path))
            submit.write("cp CONTCAR POSCAR                    \n")
            submit.write("done                                 \n")
        return jobname

    def fouropt(self, submit_dirpath):
        jobname = "fouropt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("for i in {1..4}; do                  \n")
            submit.write("cp INCAR_$i INCAR                    \n")
            submit.write("killall -9 vasp_std                  \n")
            submit.write("sleep 3                              \n")
            submit.write("{} {} > vasp.log_$i 2>&1             \n".format(self.vasp_inputpara.execmd, vaspstd_path))
            submit.write("cp CONTCAR POSCAR                    \n")
            submit.write("done                                 \n")
        return jobname

    def dfpt(self, submit_dirpath):
        jobname = "dfpt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('echo "run fine DFPT"                \n')
            submit.write('cp -f SPOSCAR POSCAR                \n')
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                        
        return jobname

    def disp(self, submit_dirpath):
        jobname = "disp.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        print(submit_script_filepath)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('echo "run Displacement" && pwd      \n')
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))   
        return jobname

    def only_onemode(self, submit_dirpath):
        jobname = "only_onemode.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))         
        return jobname
    
    def scf_eband_eledos(self, submit_dirpath):
        jobname = "scf_eband_eledos.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("cd scf\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                                                                         
            submit.write("cd ../eband\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                                                                         
            submit.write("cd ../eledos\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))                                                                         
        return jobname
    
    def scf_eband(self, submit_dirpath):
        jobname = "scf_eband.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("cd scf\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            submit.write("cd ../eband\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))       
        return jobname
    
    def scf_eledos(self, submit_dirpath):
        jobname = "scf_eledos.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("cd scf\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            submit.write("cd ../eledos\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))       
        return jobname
    
    def eband_eledos(self, submit_dirpath):
        jobname = "eband_eledos.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            submit.write("cd eband\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            submit.write("cd ../eledos\n")
            submit.write("cp ../scf/CHGCAR .\n")
            submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))       
        return jobname

    def md_nve(self, submit_dirpath):
        jobname = "md_nve.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            if vasp_inputpara.kspacing is not None:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            else:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspgam_path))                                                                       
        return jobname

    def md_nvt(self, submit_dirpath):
        jobname = "md_nvt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            if vasp_inputpara.kspacing is not None:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            else:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspgam_path))                                                                       
        return jobname
    
    def md_npt(self, submit_dirpath):
        jobname = "md_npt.sh"
        submit_script_filepath = os.path.join(submit_dirpath, jobname)
        with open(submit_script_filepath, "w") as submit:
            submit.writelines(self.jobtitle)
            if vasp_inputpara.kspacing is not None:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspstd_path))  
            else:
                submit.write('{} {} > vasp.log 2>&1               \n'.format(self.vasp_inputpara.execmd, vaspgam_path))                                                      
        return jobname