import os
import re
import logging
from pathlib import Path
from qe.qe_inputpara import qe_inputpara
from qe.qebin import qebin_path, eliashberg_x_path, bashtitle, slurmtitle, pbstitle, lsftitle


logger = logging.getLogger(__name__)

class qe_writesubmit:

    def __init__(
        self,
        qe_inputpara: qe_inputpara
        ):
 
        self.qe_inputpara = qe_inputpara

        if self.qe_inputpara.submit_job_system == "slurm":
            self.jobtitle = self.update_slurmPartition(slurmtitle, qe_inputpara.queue)
        elif self.qe_inputpara.submit_job_system == "pbs":
            self.jobtitle = pbstitle
        elif self.qe_inputpara.submit_job_system == "lsf":
            self.jobtitle = lsftitle
        elif self.qe_inputpara.submit_job_system == "bash":
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

        info = '''You can use the modes, please carefully compare your input mode is one of the above modes
        {:<20} {:<20} {:<20} {:<20}
        {:<20} {:<20} {:<20} {:<20}
        {:<20} {:<20} {:<20} {:<20}
        {:<20} {:<20}
        '''.format("relax-vc", "scffit", "scf", "prepare", "nscf", "nosplit", "split_dyn0", "split_assignQ", "q2r", "matdyn", "eletdos", "phonodos", "nscf", "McAD", "eliashberg")
        logger.debug(info)


        if mode==None:
            mode=self.qe_inputpara.mode

        if mode =="relax-vc":
            jobname = self.s1_relax(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="scffit":
            jobname = self.s2_scffit(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="scf":
            jobname = self.s3_scf(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="prepareall":
            jobname = self.s123_prepare(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="preparescf":
            jobname = self.s23_prepare(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="nosplit":
            jobname = self.s4_PhNoSplit(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="split_dyn0":
            jobnames = []
            for i, inname in enumerate(inpufilename):
                split_ph_dir = os.path.join(self.qe_inputpara.work_path, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                jobname = self.s5_PhSplitDyn0(split_ph_dir, inname)
                jobnames.append(jobname)
                logger.debug(f"finish writing submit job script in {i+1}")
            return jobnames
        if mode =="split_assignQ":
            jobnames = []
            for i, inname in enumerate(inpufilename):
                split_ph_dir = os.path.join(self.qe_inputpara.work_path, str(i+1))
                if not os.path.exists(split_ph_dir):
                    raise FileExistsError (f"There is no {split_ph_dir}")
                jobname = self.s5_PhSplitAssignQ(split_ph_dir, inname)
                jobnames.append(jobname)
                logger.debug(f"finish writing submit job script in {i+1}")
            return jobnames
        
        if mode =="q2r":
            jobname = self.s6_q2r(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="matdyn":
            jobname = self.s7_matdyn(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode == "phonobanddata":
            print(inpufilename)
            jobname = self.s7_phonobanddata(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="phonodos":
            jobname = self.s8_phonodos(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="processphono":
            jobname = self.s678_processphono(self.qe_inputpara.work_path, inpufilename)
            return jobname
        
        if mode =="eledos":
            jobname = self.s8_eledos(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="eleband":
            jobname = self.s8_eleband(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="eleproperties":
            jobname = self.s11_eleproperties(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="McAD":
            jobname = self.s9_lambda(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="eliashberg":
            jobname = self.s9_eliashberg(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="nscf":
            jobname = self.s10_nscf(self.qe_inputpara.work_path, inpufilename)
            return jobname
        
        if mode =="twin":
            jobname = self.a1_twin(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="kel":
            jobname = self.a2_kel(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="lambda_mu_k":
            jobname = self.a3_lambda_mu_k(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="scdft_tc":
            jobname = self.a4_scdft_tc(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="deltaf":
            jobname = self.a5_deltaf(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="qpdos":
            jobname = self.a6_qpdos(self.qe_inputpara.work_path, inpufilename)
            return jobname
        if mode =="sctk_all":
            jobname = self.a7_sctk_all(self.qe_inputpara.work_path, inpufilename)
            return jobname
        


    #  job scripts
    def s1_relax(self, _dirpath, inpufilename):
        _inpufilename = inpufilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s1_relax.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/pw.x -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inpufilename, _outputfilename))
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
            j.write('{} {}/pw.x -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))                                                        
        return jobname
        
    def s3_scf(self, _dirpath, inpufilename):
        _inpufilename = inpufilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s3_scf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename,  _outputfilename)) 
        return jobname

    def s123_prepare(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s123_prepare.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            relax_in, relax_out = input_output.__next__()
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, relax_in,  relax_out))
            j.write("wait\n")
            j.write("Nat=$(grep 'number of k points' -B 2 relax.out |head -n 1|awk {'print($1)'})\n")
            j.write("StruLine=$(expr $Nat + 5)\n")
            j.write("grep 'CELL_' -A $StruLine relax.out |tail -n `expr $StruLine + 1` > new_structure.out\n")
            j.write("sed  -i '/^$/d' new_structure.out\n")
            j.write("\n")
            j.write("\n")
            j.write("\n")
            j.write("\n")
            j.write("cell_parameters=$(grep -n 'CELL_PARAMETERS' scffit.in | awk -F : '{print $1}')\n")
            j.write("k_points=$(grep -n 'K_POINTS' scffit.in | awk -F : '{print $1}')\n")
            j.write("insert_position=$(expr $cell_parameters - 1)\n")
            j.write("stop_delete_position=$(expr $k_points - 1)\n")
            j.write('sed -i "${cell_parameters}, ${stop_delete_position}d" scffit.in\n')
            j.write('sed -i "${insert_position}r new_structure.out" scffit.in\n')
            scffit_in, scffit_out = input_output.__next__()
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, scffit_in,  scffit_out))
            j.write("\n")
            j.write("\n")
            j.write("\n")
            j.write("\n")
            j.write("cell_parameters=$(grep -n 'CELL_PARAMETERS' scf.in | awk -F : '{print $1}')\n")
            j.write("k_points=$(grep -n 'K_POINTS' scf.in | awk -F : '{print $1}')\n")
            j.write("insert_position=$(expr $cell_parameters - 1)\n")
            j.write("stop_delete_position=$(expr $k_points - 1)\n")
            j.write('sed -i "${cell_parameters}, ${stop_delete_position}d" scf.in\n')
            j.write('sed -i "${insert_position}r new_structure.out" scf.in\n')
            scf_in, scf_out = input_output.__next__()
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, scf_in,  scf_out))
        return jobname

    def s23_prepare(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s23_prepare.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            scffit_in, scffit_out = input_output.__next__()
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, scffit_in, scffit_out))
            scf_in, scf_out = input_output.__next__()
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, scf_in,  scf_out))
        return jobname

    def s4_PhNoSplit(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s4_PhNoSplit.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/ph.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inpufilename, _outputfilename))
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
            j.write('killall -9 vasp_std; killall -9 pw.x; killall -9 ph.x                  \n')
            j.write('echo "run scf.fit"                                                     \n')
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputscffit_name, _outputscffit_name))
            j.write('echo "run scf"                                                         \n')
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputscf_name,     _outputscf_name))
            j.write('echo "run split_ph"                                                    \n')
            j.write('{} {}/ph.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputsplitph_name,  _outputsplitph_name))   
        return jobname

    def s5_PhSplitAssignQ(self, _dirpath, inputfilename):
        _inputscffit_name, _inputscf_name, _inputsplitph_name = inputfilename
        _outputscffit_name  = _inputscffit_name.split(".")[0] + ".out"
        _outputscf_name     = _inputscf_name.split(".")[0] + ".out"
        _outputsplitph_name = _inputsplitph_name.split(".")[0] + ".out"
        jobname = "s5_PhAssignQ.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('killall -9 vasp_std; killall -9 pw.x; killall -9 ph.x                  \n')
            j.write('echo "run scf.fit"                                                \n')
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputscffit_name, _outputscffit_name))
            j.write('echo "run scf"                                                    \n')
            j.write('{} {}/pw.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputscf_name, _outputscf_name))
            j.write('echo "run split_ph"                                                    \n')
            j.write('{} {}/ph.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path,  self.qe_inputpara.npool, _inputsplitph_name,  _outputsplitph_name))   
        return jobname

    def s6_q2r(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s6_q2r.sh"
        _script_filepath = os.path.join(_dirpath,jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/q2r.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
            j.write('grep nqs q2r.out > nqs                   \n')  
        return jobname

    def s7_phonobanddata(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s7_phonobanddata.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{}/plotband.x <{}> {} \n'.format(qebin_path, _inpufilename, _outputfilename))
        return jobname
    
    def s7_matdyn(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s7_matdyn.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/matdyn.x -npool {} <{}> {} \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname
    
    def s8_eledos(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s8_eledos.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)

            nscf_in, nscf_out = input_output.__next__()
            j.write('echo "nscf"\n')
            j.write('{} {}/pw.x -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, nscf_in, nscf_out))
            
            eletdos_in, eletdos_out = input_output.__next__()
            j.write('echo "eletdos"\n')
            j.write('{} {}/dos.x -pd .true. <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path,  eletdos_in, eletdos_out))
            
            elepdos_in, elepdos_out = input_output.__next__()
            j.write('echo "elepdos"\n')
            j.write('{} {}/projwfc.x -pd .true. <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path,  elepdos_in, elepdos_out))
        return jobname

    def s8_eleband(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s8_eband.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.writelines('mkdir -p tmp/{}.save/    \n'.format(self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.charge_density_dat).absolute(),   self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.data_file_schema_xml).absolute(), self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.paw_txt).absolute(),              self.qe_inputpara.system_name))
            eleband_in, eleband_out = input_output.__next__()
            j.write('echo "eleband"\n')
            j.write('{} {}/pw.x -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, eleband_in, eleband_out))
            
            elebanddata_in, elebanddata_out = input_output.__next__()
            j.write('echo "elebanddata"\n')
            j.write('{} {}/bands.x -pd .true. <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, elebanddata_in, elebanddata_out))
            
            elebandprojdata_in, elebandprojdata_out = input_output.__next__()
            j.write('echo "elebandprojdata"\n')
            j.write('{} {}/projwfc.x -pd .true. <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, elebandprojdata_in, elebandprojdata_out))
        return jobname

    def s8_phonodos(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s8_phonodos.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/matdyn.x  <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, _inpufilename, _outputfilename))
        return jobname

    def s678_processphono(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s678_processphono.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.writelines(self.jobtitle)
            q2r_in, q2r_out = input_output.__next__()
            j.write('echo "run q2r"                   \n')
            j.write('{} {}/q2r.x -npool {} <{}> {}    \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, q2r_in, q2r_out))
            j.write('grep nqs q2r.out > nqs           \n')
            matdyn_in, matdyn_out = input_output.__next__()
            j.write('echo "run phonoband"             \n')
            j.write('{} {}/matdyn.x  <{}> {}          \n'.format(self.qe_inputpara.execmd, qebin_path, matdyn_in, matdyn_out))
            j.write('mkdir -p phonoband               \n')
            j.write('cp {}.freq     phonoband         \n'.format(self.qe_inputpara.system_name))
            j.write('cp {}.freq.gp  phonoband         \n'.format(self.qe_inputpara.system_name))
            j.write('cp {}.dyn0     phonoband         \n'.format(self.qe_inputpara.system_name))
            j.write('cp gam.lines   phonoband         \n')
            j.write('cp elph.gamma* phonoband         \n')
            plotband_in, plotband_out = input_output.__next__()
            j.write('echo "run phonoband-data-process"\n')
            j.write('{}/plotband.x <{}> {} \n'.format(qebin_path, plotband_in, plotband_out))
            phonodos_in, phonodos_out = input_output.__next__()
            j.write('echo "run phonodos"              \n')
            j.write('{} {}/matdyn.x  <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, phonodos_in, phonodos_out))
            j.write('mkdir -p phonodos                \n')
            j.write('cp {}"_phono.dos"  phonodos      \n'.format(self.qe_inputpara.system_name))
        return jobname
    
    def s9_lambda(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s9_lambda.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/lambda.x  <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, _inpufilename, _outputfilename))
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
            j.writelines('mkdir -p tmp/{}.save/    \n'.format(self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.charge_density_dat).absolute(), self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.data_file_schema_xml).absolute(), self.qe_inputpara.system_name))
            j.writelines("cp {}    tmp/{}.save/    \n".format(Path(self.qe_inputpara.paw_txt).absolute(), self.qe_inputpara.system_name))
            j.write('{} {}/pw.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def a1_twin(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "a1_twin.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/pw.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def a2_kel(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "a2_kel.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/sctk.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def a3_lambda_mu_k(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "a3_lambda_mu_k.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/sctk.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def a4_scdft_tc(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "s14_scdft_tc.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/sctk.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname

    def a5_deltaf(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "a5_deltaf.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/sctk.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname
    
    def a6_qpdos(self, _dirpath, inputfilename):
        _inpufilename = inputfilename
        _outputfilename = _inpufilename.split(".")[0] + ".out"
        jobname = "a6_qpdos.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            j.write('{} {}/sctk.x  -npool {} <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, _inpufilename, _outputfilename))
        return jobname
    
    def a7_sctk_all(self, _dirpath, inputfilenames):
        _inpufilename   = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "a7_sctk_all.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            
            nscf_in, nscf_out = input_output.__next__()
            j.write('echo "nscf"                      \n')
            j.write('{} {}/pw.x   -npool {} <{}>  {}  \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, nscf_in, nscf_out))

            twin_in, twin_out = input_output.__next__()
            j.write('echo "twin"                      \n')
            j.write('{} {}/pw.x   -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, twin_in, twin_out))
            
            kel_in, kel_out = input_output.__next__()
            j.write('echo "kel"                       \n')
            j.write('{} {}/sctk.x -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, kel_in, kel_out))
            
            lambda_mu_k_in, lambda_mu_k_out = input_output.__next__()
            j.write('echo "lambda_mu_k"               \n')
            j.write('{} {}/sctk.x -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, lambda_mu_k_in, lambda_mu_k_out))
            
            scdft_tc_in, scdft_tc_out = input_output.__next__()
            j.write('echo "scdft_tc"                  \n')
            j.write('{} {}/sctk.x -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, scdft_tc_in, scdft_tc_out))
            
            deltaf_in, deltaf_out = input_output.__next__()
            j.write('echo "deltaf"                    \n')
            j.write('{} {}/sctk.x -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, deltaf_in, deltaf_out))
            
            qpdos_in, qpdos_out = input_output.__next__()
            j.write('echo "qpdos"                     \n')
            j.write('{} {}/sctk.x -npool {} <{}> {}   \n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, qpdos_in, qpdos_out))
        return jobname
    
    def s11_eleproperties(self, _dirpath, inputfilenames):
        _inpufilename = inputfilenames
        _outputfilename = [ipname.split(".")[0] + ".out" for ipname in _inpufilename]
        input_output    = zip(_inpufilename, _outputfilename)
        jobname = "s11_eleproperties.sh"
        _script_filepath = os.path.join(_dirpath, jobname)
        with open(_script_filepath, "w") as j:
            j.writelines(self.jobtitle)
            # 先做能带计算 输入文件是 eleband.in
            eleband_in, eleband_out = input_output.__next__()
            j.write('echo "eleband"\n')
            j.write('{} {}/pw.x <{}> {}  \n\n'.format(self.qe_inputpara.execmd, qebin_path, self.qe_inputpara.npool, eleband_in, eleband_out))
            
            # 先做处理能带数据 输入文件是 elebanddata.in
            elebanddata_in, elebanddata_out = input_output.__next__()
            j.write('echo "elebanddata"\n')
            j.write('{} {}/bands.x -pd .true. <{}> {}  \n\n'.format(self.qe_inputpara.execmd, qebin_path, elebanddata_in, elebanddata_out))
            
            elebandprojdata_in, elebandprojdata_out = input_output.__next__()
            j.write('echo "elebandprojdata"\n')
            j.write('{} {}/projwfc.x      <{}> {}  \n'.format(self.qe_inputpara.execmd, qebin_path, elebandprojdata_in, elebandprojdata_out))

            # 再做非自洽计算
            nscf_in, nscf_out = input_output.__next__()
            j.write('echo "nscf"\n')
            j.write('{} {}/pw.x -npool {} <{}> {} \n\n'.format(self.qe_inputpara.execmd, qebin_path, nscf_in, nscf_out))
        
            # 再做总dos计算TDOS
            eletdos_in, eletdos_out = input_output.__next__()
            j.write('echo "eletdos"\n')
            # -pd .true. 避免调用mpirun时出错，似乎 dos.x并不能并行调用mpirun
            j.write('{} {}/dos.x -pd .true. <{}> {} \n\n'.format(self.qe_inputpara.execmd, qebin_path, eletdos_in, eletdos_out))

            # 再做投影dos计算PDOS
            elepdos_in, elepdos_out = input_output.__next__()
            j.write('echo "elepdos"\n')
            # -pd .true. 避免调用mpirun时出错，似乎 projwfc.x并不能并行调用mpirun
            j.write('{} {}/projwfc.x -pd .true. <{}> {}  \n\n'.format(self.qe_inputpara.execmd, qebin_path, elepdos_in, elepdos_out))

        return jobname