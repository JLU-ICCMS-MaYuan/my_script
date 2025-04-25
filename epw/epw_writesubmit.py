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
            self.jobtitle = self.update_slurmPartition(slurmtitle, epw_inputpara.queue)
        elif self.epw_inputpara.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.epw_inputpara.submit_job_system == "lsf":
            self.jobtitle = lsftitle
        elif self.epw_inputpara.submit_job_system == "bash":
            self.jobtitle = bashtitle
        else:
            self.jobtitle = ''
            
    def update_slurmPartition(self, title:str, new_partition:str):
        if self.check_partition_exists(new_partition):
            # 使用正则表达式替换 --partition 后面的内容
            updated_title = re.sub(r'--partition=\S+', f'--partition={new_partition}', title)
            logger.info(f"Partition exist! {new_partition}")
            return updated_title
        else:
            logger.info(f"{new_partition} doesn't exist! Keep partition name in ~/.my_scripts.py")
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
        if mode == "epw_elph":
            jobname = self.j3_epw_elph(self.epw_inputpara.work_path, inpufilename)
            return jobname
        if mode == "epw_sc":
            jobname = self.j4_epw_aniso_sc(self.epw_inputpara.work_path, inpufilename)
            return jobname

 
    def j1_epw_energyband(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j1_epw_energyband.sh"
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
            j.write('echo "epw_phono"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_phono_in, epw_phono_out))
            
            epw_phonodata_in, epw_phonodata_out = next(input_output)
            j.write('echo "epw_phonodata"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_phonodata_in, epw_phonodata_out))
            
            epw_phonodata_plot_in, epw_phonodata_plot_out = next(input_output)
            j.write('echo "epw_phonondat_plot"\n')
            j.write('{}/plotband.x <{}   \n'.format(qebin_path, epw_phonodata_plot_in))
            j.write('echo "You can plot pictures with freq.dat(phonon) and band.dat(eletron)"')
        return jobname

    def j3_epw_elph(self,  _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "j3_epw_elph.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def j4_epw_aniso_sc(self,  _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "j4_epw_sc.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            epw_iso_sc_in, epw_iso_sc_out = next(input_output)
            j.write('echo "epw_iso_sc"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_iso_sc_in, epw_iso_sc_out))
        
            epw_aniso_sc_in, epw_aniso_sc_out = next(input_output)
            j.write('echo "epw_aniso_sc"\n')
            j.write('{} {}/epw.x -npool {} <{}> {}  \n'.format(self.epw_inputpara.execmd, epwbin_path, self.epw_inputpara.npool, epw_aniso_sc_in, epw_aniso_sc_out))

        return jobname