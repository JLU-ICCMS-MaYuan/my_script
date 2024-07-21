import os
import re
import sys
import logging
import shutil
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain


from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor


from vasp.config import config
from vasp.vasp_base import vasp_base
from vasp.vasp_inputpara import *
from vasp.vasp_writeincar import vasp_writeincar
from vasp.vasp_writesubmit import vasp_writesubmit
from vasp.vasp_submitjob import vasp_submitjob

logger = logging.getLogger(__name__)


class vasp_relax:

    def __init__(self, args: ArgumentParser):
        # read input para
        _config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.relax_inputpara  = vasp_inputpara.init_from_config1(_config)

        # init the INCAR
        self._vasp_writeincar = vasp_writeincar.init_from_relaxinput(self.relax_inputpara)
        self._vasp_writeincar.writeinput()
        
        # init the submit job script
        _vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(self.relax_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        # submit the job
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(self.relax_inputpara)
        if self.relax_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)


class vaspbatch_relax:
    def __init__(self, args: ArgumentParser):
        _config = config(args).read_config()
        input_dir_path = Path(_config['input_file_path'])
        if input_dir_path.is_dir():
            input_files_path  = [input_dir_path.joinpath(pathname) for pathname in os.listdir(input_dir_path) if input_dir_path.joinpath(pathname).suffix == ".cif" or input_dir_path.joinpath(pathname).suffix == ".vasp"]
            if not input_files_path:
                print("Note: --------------------")
                print(f"    The program didn't get any structures from {input_dir_path}")
                sys.exit(1)
            work_path         = _config['work_path']        ; del _config['work_path']
            press             = _config['press']            ; del _config['press']
            submit_job_system = _config['submit_job_system']; del _config['submit_job_system']
            pass                                            ; del _config['input_file_path']
            mode              = _config['mode']             ; del _config['mode']
            for input_file_path in input_files_path:
                # prepare the POSCAR POTCAR  
                print("\nNote: --------------------")
                print(f"    Create directory for {input_file_path} file !!!")
                self.relax_inputpara  = vaspbatch_inputpara(
                    work_path=work_path,
                    press=press,
                    submit_job_system=submit_job_system,
                    input_file_path=input_file_path,
                    mode=mode,
                    **_config
                    )

                # init the INCAR
                self._vasp_writeincar  = vasp_writeincar.init_from_relaxinput(self.relax_inputpara)
                self._vasp_writeincar.writeinput()
                # init the submit job script
                _vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(self.relax_inputpara)
                jobname = _vasp_writesubmit.write_submit_scripts()
                # submit the job
                _vasp_submitjob   = vasp_submitjob.init_from_relaxinput(self.relax_inputpara)
                if self.relax_inputpara.queue is not None:
                    _vasp_submitjob.submit_mode1(jobname)
        else:
            print("Note: --------------------")
            print(f"    The {input_dir_path} path doesn't exist!")
            sys.exit(1)


class vasp_phono:
    def __init__(self, args: ArgumentParser):

        # read input para
        _config = config(args).read_config()
        
        # prepare the POSCAR POTCAR  
        phono_inputpara  = vasp_phonopara.init_from_config1(_config)

        # init the INCAR
        self._vasp_writeincar  = vasp_writeincar.init_from_phonoinput(phono_inputpara)
        self._vasp_writeincar.writeinput()

        # init the KPOINTS
        if phono_inputpara.kdensity is not None:
            phono_inputpara.create_kpoints_by_pymatgen(
                phono_inputpara.sposcar_struct_type,
                phono_inputpara.work_path.joinpath("KPOINTS"),
                phono_inputpara.kdensity,
                )
        elif phono_inputpara.kspacing is not None:
            supercell_lattice = phono_inputpara.sposcar_struct_type.lattice.matrix
            phono_inputpara.write_evenly_kpoints(
                lattice = supercell_lattice,
                kspacing=phono_inputpara.kspacing, 
                kpoints_path=phono_inputpara.work_path,
                )


        # init the submit job script
        _vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(phono_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        # submit the job
        _vasp_submitjob   = vasp_submitjob.init_from_phonoinput(phono_inputpara)
        if phono_inputpara.queue is not None:
            if phono_inputpara.mode == 'disp':
                _vasp_submitjob.submit_mode2(jobname)
            elif  phono_inputpara.mode == 'dfpt':
                _vasp_submitjob.submit_mode1(jobname)


class vaspbatch_phono:
    def __init__(self, args: ArgumentParser):

        _config = config(args).read_config()
        input_dir_path = Path(_config['input_file_path'])
        if input_dir_path.is_dir():
            input_files_path = list(input_dir_path.glob("*.vasp"))
            work_path         = _config['work_path']        ; del _config['work_path']
            press             = _config['press']            ; del _config['press']
            submit_job_system = _config['submit_job_system']; del _config['submit_job_system']
            pass                                            ; del _config['input_file_path']
            mode              = _config['mode']             ; del _config['mode']
            for input_file_path in input_files_path:
                # prepare the POSCAR POTCAR  
                phono_inputpara  = vaspbatch_phonopara(
                    work_path=work_path,
                    press=press,
                    submit_job_system=submit_job_system,
                    input_file_path=input_file_path,
                    mode=mode,
                    **_config
                    )

                # init the INCAR
                self._vasp_writeincar  = vasp_writeincar.init_from_phonoinput(phono_inputpara)
                self._vasp_writeincar.writeinput()
                # init the KPOINTS
                phono_inputpara.create_kpoints_by_pymatgen(
                    phono_inputpara.sposcar_struct_type,
                    phono_inputpara.work_path.joinpath("KPOINTS"),
                    phono_inputpara.kdensity,
                    )
                # init the submit job script
                _vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(phono_inputpara)
                jobname = _vasp_writesubmit.write_submit_scripts()
                # submit the job
                _vasp_submitjob   = vasp_submitjob.init_from_phonoinput(phono_inputpara)
                if phono_inputpara.queue is not None:
                    _vasp_submitjob.submit_mode2(jobname)


class vasp_eletron:

    def __init__(self, args: ArgumentParser):
        
        # read input para
        _config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.eletron_inputpara = vasp_eletronpara.init_from_config1(_config)

        # 准备输入文件
        self._vasp_writeincar  = vasp_writeincar.init_from_eletron(self.eletron_inputpara)
        if 'scf'    in self.eletron_inputpara.mode:
            # 准备输入文件
            scf_path = self.scf(self.eletron_inputpara.kspacing)
        if 'eband'  in self.eletron_inputpara.mode:
            # 准备输入文件
            eband_path = self.eband()
        if 'eledos' in self.eletron_inputpara.mode:
            # 准备输入文件
            eledos_path = self.eledos(self.eletron_inputpara.kspacing/2)
        if 'cohp'   in self.eletron_inputpara.mode:
            # 准备输入文件
            cohp_path = self.cohp(self.eletron_inputpara.kspacing)

        
        # 准备提交任务的脚本
        ###########################同时进行计算的任务, 作业脚本放在work_path##############################
        # 同时进行scf, eledos, eband计算
        if  ('scf'         in self.eletron_inputpara.mode) and \
            ('eledos'      in self.eletron_inputpara.mode) and \
            ('eband'       in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="scf-eband-eledos")
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname)
        # 只进行scf, eband计算
        elif ('scf'        in self.eletron_inputpara.mode) and \
             ('eledos' not in self.eletron_inputpara.mode) and \
             ('eband'      in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(mode="scf-eband")
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname)
        # 只进行scf, eledos计算
        elif ('scf'        in self.eletron_inputpara.mode) and \
             ('eledos'     in self.eletron_inputpara.mode) and \
             ('eband'  not in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(mode="scf-eledos")
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname)
        # 只进行eband, eledos计算
        elif ('scf'    not in self.eletron_inputpara.mode) and \
             ('eledos'     in self.eletron_inputpara.mode) and \
             ('eband'      in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="eband-eledos")
            chgcar_src = self.eletron_inputpara.work_path.joinpath("scf", "CHGCAR")
            if not chgcar_src.exists():
                print(f"The CHGCAR is not found in path \n{chgcar_src.absolute()}")
                print("So The program will exit")
                sys.exit(1)
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname)
        ###########################################################################
        
        
        ##################单独进行计算的任务, 作业脚本放在work_path/sub_workpath##############################
        # 只进行scf计算
        elif ('scf'        in self.eletron_inputpara.mode) and \
             ('eledos' not in self.eletron_inputpara.mode) and \
             ('eband'  not in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="only-scf",
                submitjob_path=scf_path)
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname, submit_path=scf_path)
        # 只进行eband计算
        elif ('scf'    not in self.eletron_inputpara.mode) and \
             ('eledos' not in self.eletron_inputpara.mode) and \
             ('eband'      in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="only-eband",
                submitjob_path=eband_path)
            chgcar_src = self.eletron_inputpara.work_path.joinpath("scf", "CHGCAR")
            chgcar_dst = self.eletron_inputpara.work_path.joinpath("eband", "CHGCAR")
            if not chgcar_src.exists():
                print(f"The CHGCAR is not found in path \n{chgcar_src.absolute()}")
                print("So The program will exit")
                sys.exit(1)
            else:
                shutil.copy(chgcar_src, chgcar_dst)
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname, submit_path=eband_path)
        # 只进行eledos计算
        elif ('scf'    not in self.eletron_inputpara.mode) and \
             ('eledos'     in self.eletron_inputpara.mode) and \
             ('eband'  not in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="only-eledos",
                submitjob_path=eledos_path)
            chgcar_src = self.eletron_inputpara.work_path.joinpath("scf", "CHGCAR")
            chgcar_dst = self.eletron_inputpara.work_path.joinpath("eledos", "CHGCAR")
            if not chgcar_src.exists():
                print(f"The CHGCAR is not found in path \n{chgcar_src.absolute()}")
                print("So The program will exit")
                sys.exit(1)
            else:
                shutil.copy(chgcar_src, chgcar_dst)
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname, submit_path=eledos_path)
        # 只进行cohp计算
        
        if   ('cohp'       in self.eletron_inputpara.mode):
            # 准备任务脚本
            _vasp_writesubmit = vasp_writesubmit.init_from_eletron(self.eletron_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts(
                mode="only-cohp",
                submitjob_path=cohp_path
                )
            # 提交任务
            _vasp_submitjob = vasp_submitjob.init_from_eletroninput(self.eletron_inputpara)
            if self.eletron_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(
                    jobname,
                    submit_path=cohp_path
                    )


    def scf(self, kspacing):  
        scf_path = self.eletron_inputpara.work_path.joinpath("scf")
        if not scf_path.exists():
            os.mkdir(scf_path)
        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            scf_path,
            )
        self.eletron_inputpara.get_potcar(scf_path)
        # 准备计算电子自洽的INCAR
        self._vasp_writeincar.writeinput(mode="scf", incar_path=scf_path)
        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing, 
            scf_path,
            )
        return scf_path

    def eband(self):
        scf_path = self.eletron_inputpara.work_path.joinpath("scf")
        eband_path = self.eletron_inputpara.work_path.joinpath('eband')
        if not scf_path.exists():
            os.mkdir(scf_path)
        if not eband_path.exists():
            os.mkdir(eband_path)
        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            eband_path)
        # 准备POTCAR
        self.eletron_inputpara.get_potcar(eband_path)
        # 为计算电子band准备高对称路径
        self.eletron_inputpara.write_highsymmetry_kpoints(
            self.eletron_inputpara.ase_type, 
            kpoints_path=eband_path,
            )
        # 准备计算电子band的INCAR
        self._vasp_writeincar.writeinput(mode="eband", incar_path=eband_path)
        return eband_path
        
    def eledos(self, kspacing):
        print("NOTES:  ------------------------------ ")
        print("    KSPACING in `eledos` have to be twice than that in `scf` ")
        scf_path    = self.eletron_inputpara.work_path.joinpath("scf")
        eledos_path = self.eletron_inputpara.work_path.joinpath('eledos')
        if not scf_path.exists():
            os.mkdir(scf_path)
        if not eledos_path.exists():
            os.mkdir(eledos_path)   
        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            eledos_path)    
        self.eletron_inputpara.get_potcar(eledos_path)
        # 准备计算电子DOS的INCAR
        self._vasp_writeincar.writeinput(mode='eledos', incar_path=eledos_path)
        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing, 
            kpoints_path=eledos_path,
            )
        return eledos_path

    def cohp(self, kspacing):
        cohp_path = self.eletron_inputpara.work_path.joinpath('cohp')
        if not cohp_path.exists():
            os.mkdir(cohp_path)  
        # 准备POSCAR
        self.eletron_inputpara.get_struct_info(
            self.eletron_inputpara.struct_type,
            cohp_path)    
        self.eletron_inputpara.get_potcar(cohp_path)
        # 准备计算电子cohp的INCAR
        self._vasp_writeincar.writeinput(mode='cohp', incar_path=cohp_path)
        # 为电子自洽均匀撒点准备KPOINTS
        self.eletron_inputpara.write_evenly_kpoints(
            self.eletron_inputpara.cell_parameters,
            kspacing, 
            kpoints_path=cohp_path,
            )
        return cohp_path


class vasp_processdata(vasp_base):

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config         = config(args).read_config()
        self.input_file_path = Path(self._config['input_file_path'])
        if self._config['work_path']:
            self.work_path   = Path(self._config['work_path'])
        else:
            self.work_path   = Path.cwd()
        self.mode            = self._config['mode']
        self.ase_type        = read(self.input_file_path)
        self.struct_type     = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        
        if self.mode == "dispband" or self.mode == "dfptband":
            self.post_progress_phono_band()
        if self.mode == "dispphdos" or self.mode == "dfptphdos":
            self.post_progress_phono_dos()
        if self.mode == "eband":
            self.post_progress_eletron_band()
        if self.mode == "eledos":
            self.post_progress_eletron_dos()
        if self.mode == "hspp":
            self.read_hspp(self.work_path.joinpath("KPOINTS"))
    
    # 绘制 phonoband 
    def post_progress_phono_band(self):

        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            print("WARNING: ----------------------------- ")
            raise ValueError("    you have to specify the supercell=[?,?,?]. If you didn't specify it, maybe somthing wrong will occur !")

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [20, 20, 20]
            print("NOTES: ------------------------------ ")
            print("    you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if self.mode == "dispband":
            print("NOTES: ------------------------------ ")
            print("    Run disp-band-post-progress-module")
            _disp_num = len(list(Path(self.work_path).glob("disp-*")))
            disp_num  = str(_disp_num).rjust(3, '0')
            cwd = os.getcwd()
            os.chdir(self.work_path)
            print("NOTES: ------------------------------ ")
            print("    Run the order `phonopy -f disp-{001..%s}/vasprun.xml`" %(disp_num))
            print("    Please confirm  disp-{001..%s} are all correctively computed !!!" %(disp_num))
            os.system("phonopy -f disp-{001..%s}/vasprun.xml" %(disp_num))
            os.chdir(cwd) 


            print("\nNote: --------------------------------")
            vaspkitflag = input("If you have installed vaspkit and you want to use it, input: Yes\n")
            if vaspkitflag:
                cwd = os.getcwd()
                os.chdir(self.work_path)
                os.system('echo -e "3\n305\n3" | vaspkit')
                shutil.copy("KPATH.phonopy", "band.conf")
                diminfo = "DIM={}".format(' '.join(list(map(str, self.supercell))))
                os.system("sed -i '2s/.*/{}/' band.conf".format(diminfo))
                os.system("sed -i '/FORCE_CONSTANTS = READ/d' band.conf")
                mpinfo = "MP={}".format(' '.join(list(map(str, self.mp))))
                os.system("sed -i '6s/.*/{}/' band.conf".format(mpinfo))
                os.chdir(cwd)
            else: 
                path_name_list, path_coords = self.get_hspp(self.ase_type)
                self.write_disp_band_conf(
                    self.work_path, 
                    self.species, 
                    self.supercell, 
                    path_name_list,
                    path_coords
                    )
                
            cwd = os.getcwd()
            os.chdir(self.work_path)
            # traditional method
            print("NOTES: ------------------------------ ")
            print("    Run the order `phonopy -p -s band.conf -c POSCAR-init`")
            print("    Please confirm the POSCAR-init exist !!!")
            os.system("phonopy -p -s band.conf -c POSCAR-init")
            print("NOTES: ------------------------------ ")
            print("    Run the order `phonopy-bandplot  --gnuplot> band.dat`")
            print("    band.dat is the data for plot phononband in Origin")
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)

        elif self.mode == "dfptband":
            print("NOTES: ------------------------------ ")
            print("    Run dfpt-band-post-progress-module")
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy --fc vasprun.xml")
            os.chdir(cwd)


            print("\nNote: --------------------------------")
            vaspkitflag = input("If you have installed vaspkit and you want to use it, input: Yes\n")
            if vaspkitflag:
                cwd = os.getcwd()
                os.chdir(self.work_path)
                os.system('echo -e "3\n305" | vaspkit')
                shutil.copy("KPATH.phonopy", "band.conf")
                diminfo = "DIM={}".format(' '.join(list(map(str, self.supercell))))
                os.system("sed -i '2s/.*/{}/' band.conf".format(diminfo))
                os.chdir(cwd)
            else: 
                special_points, path_coords = self.get_hspp(self.ase_type)
                self.write_dfpt_band_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                special_points,
                path_coords
                )     

            cwd = os.getcwd()
            os.chdir(self.work_path)
            # traditional method
            os.system("phonopy --dim='{}' -p -s band.conf -c POSCAR-init".format(' '.join(list(map(str, self.supercell)))))
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)
    
    # 绘制 phonodos 
    def post_progress_phono_dos(self):

        if "phdos" in self._config:
            self.pdos = self._config['pdos']
        else:
            self.pdos = "AUTO"

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [20, 20, 20]
            print("NOTES: ------------------------------ ")
            print("    you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            print("WARNING: ----------------------------- ")
            raise ValueError("    you have to specify the supercell=[?,?,?]. If you didn't specify it, maybe somthing wrong will occur !")

        if "tmin" in self._config:
            self.tmin =  self._config['tmin']
        else:
            self.tmin = 0

        if "tmax" in self._config:
            self.tmax =  self._config['tmax']
        else:
            self.tmax = 5000

        if "tstep" in self._config:
            self.tstep =  self._config['tstep']
        else:
            self.tstep = 100

        if self.mode == "dispphdos": 
            # 获得total_dos.dat
            self.write_disp_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
                self.tmin,
                self.tmax,   
                self.tstep,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf") # -p: dos plot   -t: thermal eletron print
            os.system("phonopy -p mesh.conf -c {}".format(self.input_file_path.name)) # 获得 total_dos.dat
            os.chdir(cwd)

            # 获得pdos.dat
            self.write_disp_phdos_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
                self.pdos,
                
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p pdos.conf -c {}".format(self.input_file_path.name)) # 获得 total_dos.dat
            os.chdir(cwd)

        elif self.mode == "dfptphdos":
            # 获得total_dos.dat
            self.write_dfpt_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
                self.tmin,
                self.tmax,   
                self.tstep,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf") # -p: dos plot   -t: thermal eletron print
            os.system("phonopy -p mesh.conf -c {}".format(self.input_file_path.name)) # 获得 total_dos.dat
            os.chdir(cwd)

            # 获得pdos.dat
            self.write_dfpt_phdos_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
                self.pdos,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p pdos.conf -c {}".format(self.input_file_path.name)) # 获得 total_dos.dat
            os.chdir(cwd)
    
    # 绘制 eband
    def post_progress_eletron_band(self):

        import matplotlib.pyplot as plt
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import BSPlotter

        # 检查费米能级
        vasprunxml_path  = self.work_path.joinpath("vasprun.xml")
        self.check_efermi_energy(vasprunxml_path)
        vasprun = Vasprun(vasprunxml_path)
        e_fermi_fromband  = vasprun.efermi
        print("NOTES: ------------------------------ ")
        print("    Check the E-fermi( {} ) is equal to e_fermi_fromscf whether or not at last time !".format(e_fermi_fromband))

        eband = vasprun.get_band_structure(line_mode=True)

        e_fermi_fromvasp = vasprun.efermi

        # set figure parameters, draw figure
        bsplotter = BSPlotter(bs=eband)
        bsplotter.bs_plot_data()
        bsplotter.get_plot()
        ebandpng_path = self.work_path.joinpath('eband.png')
        bsplotter.save_plot(
            ebandpng_path, 
            img_format='png',
            ylim=[-5,  5],
            )
    
    # 绘制 eledos
    def post_progress_eletron_dos(self):

        import matplotlib.pyplot as plt
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import DosPlotter

        # 检查费米能级
        vasprunxml_path  = self.work_path.joinpath("vasprun.xml")
        self.check_efermi_energy(vasprunxml_path)
        vasprun = Vasprun(vasprunxml_path)
        e_fermi_fromdos  = vasprun.efermi
        print("NOTES: ------------------------------ ")
        print("    Check the E-fermi( {} ) is equal to e_fermi_fromscf whether or not at last time !".format(e_fermi_fromdos))
        
        # 获得dos的数据
        try:
            eledos = vasprun.complete_dos_normalized
            print("NOTES: ------------------------------ ")
            print("    Success getting the complete_dos_normalized")        
        except:
            eledos = vasprun.complete_dos
            print("NOTES: ------------------------------ ")
            print("    Success getting the complete_dos")        


        # 处理dos的数据，并决定
        dosplotter = DosPlotter()
        self.pdostype   = self._config.get('pdostype', None)
        if self.pdostype == "ele":
            print("NOTES: ------------------------------ ")
            dosplotter.add_dos_dict(eledos.get_element_dos())
            print("    Success getting the dos projected to elements")        
        elif self.pdostype == "spd":
            print("NOTES: ------------------------------ ")
            dosplotter.add_dos_dict(eledos.get_spd_dos())
            print("    Success getting the dos projected to spd-orbits")        
        elif self.pdostype == "elespd":
            print("NOTES: ------------------------------ ")
            dosplotter.add_dos_dict(eledos.get_element_spd_dos())
            print("    Success getting the dos projected to elements and spd-orbits")        
        else:
            print("NOTES: You set nothing for pdostype --")
            dosplotter.add_dos_dict(eledos.get_element_dos())
            print("    Default value is projected to element. Success getting the dos projected to elements")        
        dosplotter.get_plot()
        eledospng_path = self.work_path.joinpath('eledos.png')
        dosplotter.save_plot(
            filename=eledospng_path, 
            img_format='png', 
            xlim=[-5,  5], 
            ylim=[0,  20],
            )
    
    # 创建band.conf  目的为了获得 band of phonon
    def write_disp_band_conf(
        self, 
        band_conf_dirpath, 
        species, 
        supercell, 
        path_name_list, 
        path_coords
        ):
        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))

        band_conf_filepath = os.path.join(band_conf_dirpath, "band.conf")
        with open(band_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("NPOINTS=101             \n")
            f.write("EIGENVECTORS=.TRUE.     \n")
            f.write("BAND_LABELS={}          \n".format(' '.join(path_name_list)))
            path_coords = list(chain.from_iterable(path_coords)); path_coords=list(map(str, path_coords))
            f.write("BAND={}                 \n".format(' '.join(path_coords)))  
    
    # 创建band.conf  目的为了获得 band of phonon
    def write_dfpt_band_conf(
        self, 
        band_conf_dirpath, 
        species, 
        supercell,  
        special_points, 
        path_coords
        ):
        
        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        
        band_conf_filepath = os.path.join(band_conf_dirpath, "band.conf")
        with open(band_conf_filepath, "w") as f:
            f.write("FORCE_CONSTANTS=READ    \n")
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("NPOINTS=101             \n")
            f.write("EIGENVECTORS=.TRUE.     \n")
            f.write("BAND_LABELS={}          \n".format(' '.join(special_points)))
            path_coords = list(chain.from_iterable(path_coords)); path_coords=list(map(str, path_coords))
            f.write("BAND={}                 \n".format(' '.join(path_coords)))
    
    # 创建mesh.conf  目的为了获得 thermal_properties.yaml
    def write_disp_mesh_conf(
        self,
        mesh_conf_dirpath, 
        species, 
        supercell,
        mp,
        tmin,
        tmax,
        tstep,
        ): 

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        mesh_conf_filepath = os.path.join(mesh_conf_dirpath, "mesh.conf")
        with open(mesh_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
            f.write("TPROP=T                 \n")
            f.write("TMIN={}                 \n".format(tmin))
            f.write("TMAX={}                 \n".format(tmax))
            f.write("TSTEP={}                \n".format(tstep))
                  
    # 创建mesh.conf  目的为了获得 thermal_properties.yaml
    def write_dfpt_mesh_conf(
        self,
        mesh_conf_dirpath, 
        species, 
        supercell,
        mp,
        tmin,
        tmax,
        tstep,
        ): 

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        mesh_conf_filepath = os.path.join(mesh_conf_dirpath, "mesh.conf")
        with open(mesh_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
            f.write("FORCE_CONSTANTS = READ  \n")
            f.write("TPROP=T                 \n")
            f.write("TMIN={}                 \n".format(tmin))
            f.write("TMAX={}                 \n".format(tmax))
            f.write("TSTEP={}                \n".format(tstep))
    
    # 创建pdos.conf  目的为了获得 pdos 
    def write_disp_phdos_conf(
        self,
        pdos_conf_dirpath, 
        species, 
        supercell,
        mp,
        pdos,
        ): 

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        pdos_conf_filepath = os.path.join(pdos_conf_dirpath, "pdos.conf")
        with open(pdos_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
            f.write("PDOS = {}               \n".format(pdos))
    
    # 创建pdos.conf  目的为了获得 pdos 
    def write_dfpt_phdos_conf(
        self,
        pdos_conf_dirpath, 
        species, 
        supercell,
        mp,
        pdos,
        ): 

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        pdos_conf_filepath = os.path.join(pdos_conf_dirpath, "pdos.conf")
        with open(pdos_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
            f.write("FORCE_CONSTANTS = READ  \n")
            f.write("PDOS = {}               \n".format(pdos))

    # 检查费米能级
    def check_efermi_energy(self, vasprunxml_path):

        scfoutcar_path   = self.work_path.absolute().parent.joinpath("scf", "OUTCAR")
        # 检查费米能级
        e_fermi_fromscf  = os.popen(f"grep E-fermi {scfoutcar_path} | tail -n 1 " + "| awk '{print $3}' ").read().strip("\n")
        e_fermi_dos_band = os.popen(f"grep efermi  {vasprunxml_path}" + "| awk '{print $3}' ").read().strip("\n")
        print("NOTES: ------------------------------ ")
        print("    You have to confirm that the Fermi energy is from scf/OUTCAR. Because the Fermi energy in dos/DOSCAR is not accurate")
        print("    You can use 'grep E-fermi scf/OUTCAR' to check the Fermi energy by yourself !")
        print("    E-fermi in scf is {}".format(e_fermi_fromscf))
        print("    E-fermi in dos is {}".format(e_fermi_dos_band))
        print("    The program will use `e_fermi_fromscf` to cover the `e_fermi_dos_band`")
        if abs(float(e_fermi_fromscf) - float(e_fermi_dos_band)) > 0.0001:
            replace_efermi_in_vasprunxml = """ sed -E -i.bak """ + \
                                           """ 's/<i name="efermi">\s*[0-9]+\.[0-9]+\s*<\/i>/<i name="efermi">    {} <\/i>/' """.format(e_fermi_fromscf) + \
                                           """ {} """.format(vasprunxml_path)
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system(replace_efermi_in_vasprunxml)
            os.chdir(cwd)
        print("NOTES: ------------------------------ ")
        print("    If you wanna plot band or dos by yourself, you'd better replace efermi in DOSCAR with that in scf/OUTCAR")


class vasp_clear:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()
        self.work_path       = Path(self._config['work_path'])
        self.mode            = self._config['mode']
        if self.mode == "all":
            self.clear_all()

    def clear_all(self):
        """
        delete all except "POSCAR", "PPOSCAR", "POTCAR", "OUTCAR", "INCAR*", "*.sh", "*.vasp", "*.slurm"
        """
        reserved_files = ["POTCAR", "OUTCAR", "CONTCAR"]
        current_files = os.listdir(self.work_path)
        for file in current_files:
            if file in reserved_files:
                pass
            elif "INCAR" in file:
                pass
            elif "POSCAR" in file:
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".sh":
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".vasp":
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".slurm":
                pass
            else:
                os.system(f"rm -fr {str(Path(self.work_path).joinpath(file))}")



