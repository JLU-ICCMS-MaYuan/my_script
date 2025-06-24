import os
import re
import logging
from pathlib import Path
from epw.epw_inputpara import epw_inputpara
from epw.epwbin import epwbin_path, epwbin_path, qebin_path, bashtitle, slurmtitle, pbstitle, lsftitle


logger = logging.getLogger(__name__)

class epw_writesubmit:

    def __init__(
        self,
        epw_inputpara: epw_inputpara
        ):
 
        self.epw_inputpara = epw_inputpara

        if self.epw_inputpara.submit_job_system == "slurm":
            self.jobtitle = self.update_slurminfo(slurmtitle, epw_inputpara.queue, epw_inputpara.nodes, epw_inputpara.ntasks, epw_inputpara.ntasks_per_node, epw_inputpara.cpus_per_task)
        elif self.epw_inputpara.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.epw_inputpara.submit_job_system == "lsf":
            self.jobtitle = self.update_lsfinfo(lsftitle, epw_inputpara.ntasks, epw_inputpara.ntasks_per_node)
        elif self.epw_inputpara.submit_job_system == "bash":
            self.jobtitle = bashtitle
        else:
            self.jobtitle = ''
            
    def update_slurminfo(self, title:str, new_partition:str, nodes:int, ntasks:int, ntasks_per_node:int, cpus_per_task:int):
        if self.check_partition_exists(new_partition):
            # 使用正则表达式替换 --partition 后面的内容
            updated_title = re.sub(r'--partition=\S+', f'--partition={new_partition}', title)
            logger.info(f"Partition exist! {new_partition}")
        else:
            logger.info(f"{new_partition} doesn't exist! Keep partition name in ~/.my_scripts.py")
            updated_title = title
            
        if nodes is not None:
            updated_title = re.sub(r'--nodes=\S+', f'--nodes={nodes}', updated_title)
        if ntasks is not None:
            updated_title = re.sub(r'--ntasks=\S+', f'--ntasks={ntasks}', updated_title)
        if ntasks_per_node is not None:
            updated_title = re.sub(r'--ntasks-per-node=\S+', f'--ntasks-per-node={ntasks_per_node}', updated_title)
        if cpus_per_task is not None:
            updated_title = re.sub(r'--cpus-per-task=\S+', f'--cpus-per-task={cpus_per_task}', updated_title)
        return updated_title
    
    def update_lsfinfo(self, title:str, ntasks:int, ntasks_per_node:int):
        updated_title = title
        if ntasks is not None:
            updated_title = re.sub(r"#BSUB -n\s+.*",    f'#BSUB -n {ntasks}', updated_title)
            # '.*' 是正则表达式中的核心：
            # .：匹配任意一个字符（除了换行符）
            # *：表示“零个或多个”前面的字符
        if ntasks_per_node is not None:
            updated_title = re.sub(r"#BSUB -R\s+'.*'" , f"#BSUB -R 'span[ptile={ntasks_per_node}]'", updated_title)
        return updated_title


    
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
            logger.warning(f"Error throws up when run `sinfo -h --format=%P`: {e}")
            return False

    def write_submit_scripts(self, inpufilename, mode=None):

        if mode==None:
            mode=self.epw_inputpara.mode
            info = 'You can use the modes, please carefully compare your input mode is one of the above modes\n' +\
            '{:<20} {:<20} {:<20} {:<20}'.format("epw_band", "epw_phono", "McAD", "eliashberg")
            logger.debug(info)
        if mode == "epw_eband":
            jobname = self.j1_epw_energyband(self.epw_inputpara.work_path, inpufilename)
            return jobname
        if mode == "epw_phono":
            jobname = self.j2_epw_phono(self.epw_inputpara.work_path, inpufilename)
            return jobname
        if mode == "epw_phonodata":
            jobname = self.j3_epw_phonodata(self.epw_inputpara.work_path, inpufilename)
            return jobname
        if mode == "epw_elph":
            jobname = self.j4_epw_elph(self.epw_inputpara.work_path, inpufilename)
            return jobname
        if mode == "epw_sc":
            if len(self.epw_inputpara.muc) > 1:
                jobname = []
                for mu in self.epw_inputpara.muc:
                    iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(mu))
                    iso_jobname = self.j5_epw_iso_sc(iso_mu_path, inpufilename[0])
                    aniso_mu_path = self.epw_inputpara.work_path.joinpath("aniso_muc_{}".format(mu))
                    aniso_jobname = self.j6_epw_aniso_sc(aniso_mu_path, inpufilename[1])
                return iso_jobname, aniso_jobname
            else:
                iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(self.epw_inputpara.muc))
                iso_jobname = self.j5_epw_iso_sc(iso_mu_path, inpufilename[0])
                aniso_mu_path = self.epw_inputpara.work_path.joinpath("aniso_muc_{}".format(self.epw_inputpara.muc))
                aniso_jobname = self.j6_epw_aniso_sc(aniso_mu_path, inpufilename[1])
                return iso_jobname, aniso_jobname
        if mode == "epw_prtgkk":
            prtgkk_path = self.epw_inputpara.work_path.joinpath("prtgkk")
            jobname = self.j7_epw_prtgkk(prtgkk_path, inpufilename)
            return jobname
        if mode == "epw_fermi_nest":
            fermi_nest = self.epw_inputpara.work_path.joinpath("fermi_nest")
            jobname = self.j8_epw_fermi_nest(fermi_nest, inpufilename)
            return jobname
        if mode == "epw_linearized_iso":
            if len(self.epw_inputpara.muc) > 1:
                jobname = []
                for mu in self.epw_inputpara.muc:
                    iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(mu))
                    iso_jobname = self.j9_epw_linearized_iso(iso_mu_path, inpufilename[0])
                return iso_jobname
            else:
                iso_mu_path = self.epw_inputpara.work_path.joinpath("iso_muc_{}".format(self.epw_inputpara.muc))
                iso_jobname = self.j9_epw_linearized_iso(iso_mu_path, inpufilename[0])
                return iso_jobname
    def j1_epw_energyband(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j1_epw_eband.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def j2_epw_phono(self,  _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "j2_epw_phono.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            
            epw_phono_in, epw_phono_out = next(input_output)
            j.write('echo "epw_phonodata"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_phono_in, epw_phono_out))
            
            epw_phonodata_plot_in, epw_phonodata_plot_out = next(input_output)
            j.write('echo "epw_phonondat_plot"\n')
            j.write('{}/plotband.x <{}   \n'.format(qebin_path, epw_phonodata_plot_in))
            j.write('echo "You can plot pictures with freq.dat(phonon) and band.dat(eletron)"')
        return jobname

    def j3_epw_phonodata(self,  _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "j3_epw_phonodata.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            epw_phonodata_in, epw_phonodata_out = next(input_output)
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_phonodata_in, epw_phonodata_out))
            
            epw_phonodata_plot_in, epw_phonodata_plot_out = next(input_output)
            j.write('echo "epw_phonondat_plot"\n')
            j.write('{}/plotband.x <{}   \n'.format(qebin_path, epw_phonodata_plot_in))
            j.write('echo "You can plot pictures with freq.dat(phonon) and band.dat(eletron)"')

        return jobname

    def j4_epw_elph(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j4_epw_elph.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def j5_epw_iso_sc(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j5_epw_iso_sc.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('cp {}.a2f backup_{}.a2f\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
            j.write('cp {}.a2f backup_{}.a2f_proj\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
            if self.epw_inputpara.keep_a2f:
                j.write('echo "romove a2f a2f_proj"\n')
                j.write('rm {}.a2f {}.a2f_proj\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
                j.write('abspath=`realpath ..`\n')
                j.write('mkdir -p ./tmp/{}.ephmat\n'.format(self.epw_inputpara.system_name))  # 保证目录存.format(self.epw)在
                j.write('ln -sf ${abspath}' + '/tmp/{}.ephmat/* ./tmp/{}.ephmat/\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}' + '/tmp/{}.epmatwp  ./tmp/\n'.format(self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}' + '/{}.ukk          . \n'.format(self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}/crystal.fmt            .\n')
                j.write('ln -sf ${abspath}/restart.fmt            .\n')
                j.write('ln -sf ${abspath}/selecq.fmt             .\n')
            j.write('echo "epw_iso_sc"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        
        return jobname
    
    def j6_epw_aniso_sc(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j6_epw_aniso_sc.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('cp {}.a2f backup_{}.a2f\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
            j.write('cp {}.a2f backup_{}.a2f_proj\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
            if self.epw_inputpara.keep_a2f:
                j.write('echo "romove a2f a2f_proj"\n')
                j.write('rm {}.a2f {}.a2f_proj\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
                j.write('abspath=`realpath ..`\n')
                j.write('mkdir -p ./tmp/{}.ephmat\n'.format(self.epw_inputpara.system_name))  # 保证目录存.format(self.epw)在
                j.write('ln -sf ${abspath}' + '/tmp/{}.ephmat/* ./tmp/{}.ephmat/\n'.format(self.epw_inputpara.system_name, self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}' + '/tmp/{}.epmatwp  ./tmp/\n'.format(self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}' + '/{}.ukk          . \n'.format(self.epw_inputpara.system_name))
                j.write('ln -sf ${abspath}/crystal.fmt            .\n')
                j.write('ln -sf ${abspath}/restart.fmt            .\n')
                j.write('ln -sf ${abspath}/selecq.fmt             .\n')
            j.write('echo "epw_aniso_sc"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))

        return jobname
    
    def j7_epw_prtgkk(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j7_epw_prtgkk.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('echo "calcualte el-ph matrix elements"\n')
            j.write('abspath=`realpath ..`\n')
            j.write('mkdir -p ./tmp\n')  
            j.write('ln -sf ${abspath}' + '/tmp/{}.epmatwp  ./tmp/\n'.format(self.epw_inputpara.system_name))
            j.write('ln -sf ${abspath}' + '/{}.ukk          . \n'.format(self.epw_inputpara.system_name))
            j.write('ln -sf ${abspath}/crystal.fmt            .\n')
            j.write('ln -sf ${abspath}/restart.fmt            .\n')
            j.write('ln -sf ${abspath}/epwdata.fmt            .\n')
            j.write('ln -sf ${abspath}/vmedata.fmt            .\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname
    
    def j8_epw_fermi_nest(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j8_epw_fermi_nest.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('echo "calcualte fermi nest"\n')
            j.write('abspath=`realpath ..`\n')
            j.write('mkdir -p ./tmp\n')  
            j.write('ln -sf ${abspath}' + '/tmp/{}.epmatwp  ./tmp/\n'.format(self.epw_inputpara.system_name))
            j.write('ln -sf ${abspath}' + '/{}.ukk          . \n'.format(self.epw_inputpara.system_name))
            j.write('ln -sf ${abspath}/crystal.fmt            .\n')
            j.write('ln -sf ${abspath}/epwdata.fmt            .\n')
            j.write('ln -sf ${abspath}/vmedata.fmt            .\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def j9_epw_linearized_iso(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j9_epw_linearized_iso.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('echo "epw_linearized_iso"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname
