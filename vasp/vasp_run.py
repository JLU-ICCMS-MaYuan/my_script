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


def vasp_relax(args: ArgumentParser) -> None:
    # read input para
    _config = config(args).read_config()

    # prepare the POSCAR POTCAR  
    relax_inputpara  = vasp_inputpara.init_from_config1(_config)

    # init the INCAR
    _vasp_writeincar  = vasp_writeincar.init_from_relaxinput(relax_inputpara)
    if relax_inputpara.mode == 'rvf' or relax_inputpara.mode == 'rv1':
        _vasp_writeincar.opt_fine_incar(relax_inputpara.work_path) 
    elif relax_inputpara.mode == 'rv3':
        _vasp_writeincar.opt_incar1(relax_inputpara.work_path)
        _vasp_writeincar.opt_incar2(relax_inputpara.work_path)
        _vasp_writeincar.opt_incar3(relax_inputpara.work_path)
    # init the submit job script
    _vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(relax_inputpara)
    jobname = _vasp_writesubmit.write_submit_scripts()
    # submit the job
    _vasp_submitjob = vasp_submitjob.init_from_relaxinput(relax_inputpara)
    if relax_inputpara.queue is not None:
        _vasp_submitjob.submit_mode1(jobname)


def vaspbatch_relax(args: ArgumentParser) -> None:

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
            print("\nNote: --------------------")
            print("    Create directory for {} file !!!")
            relax_inputpara  = vaspbatch_inputpara(
                work_path=work_path,
                press=press,
                submit_job_system=submit_job_system,
                input_file_path=input_file_path,
                mode=mode,
                **_config
                )

            # init the INCAR
            _vasp_writeincar  = vasp_writeincar.init_from_relaxinput(relax_inputpara)
            if relax_inputpara.mode == 'rvf' or relax_inputpara.mode == 'rv1':
                _vasp_writeincar.opt_fine_incar(relax_inputpara.work_path) 
            elif relax_inputpara.mode == 'rv3':
                _vasp_writeincar.opt_incar1(relax_inputpara.work_path)
                _vasp_writeincar.opt_incar2(relax_inputpara.work_path)
                _vasp_writeincar.opt_incar3(relax_inputpara.work_path)
            # init the submit job script
            _vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(relax_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts()
            # submit the job
            _vasp_submitjob   = vasp_submitjob.init_from_relaxinput(relax_inputpara)
            if relax_inputpara.queue is not None:
                _vasp_submitjob.submit_mode1(jobname)


def vasp_phono(args: ArgumentParser) -> None:

    # read input para
    _config = config(args).read_config()
    
    # prepare the POSCAR POTCAR  
    phono_inputpara  = vasp_phonopara.init_from_config1(_config)

    # init the INCAR
    _vasp_writeincar  = vasp_writeincar.init_from_phonoinput(phono_inputpara)
    if phono_inputpara.mode == 'disp':
        _vasp_writeincar.disp_incar(phono_inputpara.work_path)
    elif phono_inputpara.mode == 'dfpt':
        _vasp_writeincar.dfpt_incar(phono_inputpara.work_path)

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


def vaspbatch_phono(args: ArgumentParser) -> None:

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
            _vasp_writeincar  = vasp_writeincar.init_from_phonoinput(phono_inputpara)
            if phono_inputpara.mode == 'disp':
                _vasp_writeincar.disp_incar(phono_inputpara.work_path)
            elif phono_inputpara.mode == 'dfpt':
                _vasp_writeincar.dfpt_incar(phono_inputpara.work_path)
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


def vasp_eletron(args):
        
    # read input para
    _config = config(args).read_config()

    # prepare the POSCAR POTCAR  
    eletron_inputpara = vasp_inputpara.init_from_config1(_config)

    # init the INCAR
    _vasp_writeincar  = vasp_writeincar.init_from_eletron(eletron_inputpara)
    if eletron_inputpara.mode == 'scf':
        _vasp_writeincar.scf_incar(
            eletron_inputpara.work_path
            )
        eletron_inputpara.write_evenly_kpoints(
            eletron_inputpara.kspacing, 
            eletron_inputpara.work_path
            )
        _vasp_writesubmit = vasp_writesubmit.init_from_eletron(eletron_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(eletron_inputpara)
        if eletron_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
    elif eletron_inputpara.mode == 'eband':
        chgcar_src = eletron_inputpara.work_path.parent.joinpath("scf", "CHGCAR")
        chgcar_dst = eletron_inputpara.work_path.joinpath("CHGCAR")
        if chgcar_src.exists():
            shutil.copy(chgcar_src, chgcar_dst)
        else:
            print("NOTES: ------------------------------ ")
            print(f"    The CHGCAR doesn't exist in {chgcar_src}")
            sys.exit(1)
        eletron_inputpara.write_highsymmetry_kpoints(
            eletron_inputpara.ase_type, 
            eletron_inputpara.work_path.joinpath("KPOINTS"),
            )
        _vasp_writeincar.band_incar(eletron_inputpara.work_path)
        _vasp_writesubmit = vasp_writesubmit.init_from_eletron(eletron_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(eletron_inputpara)
        if eletron_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
        print("NOTES: ------------------------------ ")
        print("If you meet the erros in eband just like: ")
        print("    ERROR: charge density could not be read from file CHGCAR for ICHARG>10")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ in scf/OUTCAR are different from those in eband/OUTCAR ")
        print("    SOLUTION: You can let NGX,NGY,NGZ in eledos/INCAR be the same as scf/OUTCAR")
        print("If you meet the erros in eledos just like: ")
        print("    WARING: Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accurate calculation. Thus, the results might be wrong. ")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ you'v customized aren't matched with the PREC=Accurate ")
        print("    SOLUTION: You can let PREC=Normal eband/INCAR")
    elif eletron_inputpara.mode == 'eledos':
        chgcar_src = eletron_inputpara.work_path.parent.joinpath("scf", "CHGCAR")
        chgcar_dst = eletron_inputpara.work_path.joinpath("CHGCAR")
        if chgcar_src.exists():
            shutil.copy(chgcar_src, chgcar_dst)
        else:
            print("NOTES: ------------------------------ ")
            print(f"    The CHGCAR doesn't exist in {chgcar_src}")
            sys.exit(1)
        print("NOTES:  ------------------------------ ")
        print("    KSPACING in `eledos` have to be twice than that in `scf` ")
        eletron_inputpara.write_evenly_kpoints(
            eletron_inputpara.kspacing, 
            eletron_inputpara.work_path
            )
        _vasp_writeincar.eledos_incar(eletron_inputpara.work_path)
        _vasp_writesubmit = vasp_writesubmit.init_from_eletron(eletron_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(eletron_inputpara)
        if eletron_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
        print("NOTES: ------------------------------ ")
        print("if you meet the erros just like: ")
        print("    ERROR: charge density could not be read from file CHGCAR for ICHARG>10")
        print("    Possible reason is that NGX, NGY, NGZ in scf/OUTCAR are different from those in eledos/OUTCAR ")
        print("    If you wanna to solve it, you need let NGX,NGY,NGZ in eledos/INCAR be the same as scf/OUTCAR")
        print("If you meet the erros just like: ")
        print("    WARING: Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accurate calculation. Thus, the results might be wrong. ")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ you'v customized aren't matched with the PREC=Accurate ")
        print("    SOLUTION: You can let PREC=Normal in eledos/INCAR ")
    else:
        print("NOTES: ------------------------------ ")
        print(f"    The mode {vasp_writeincar.mode} you specified isn't supported. ")


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
        if self.mode == "dispphdos" or self.mode == "dispphdos":
            self.post_progress_phono_dos()
        if self.mode == "eband":
            self.post_progress_eletron_band()
        if self.mode == "eledos":
            self.post_progress_eletron_dos()
        if self.mode == "hspp":
            self.get_hspp(self.ase_type)
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
            self.mp = [8,8,8]
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

        if self.mode == "dispphdos": 
            # 获得total_dos.dat
            self.write_disp_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
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

        elif self.mode == "dispphdos":
            # 获得total_dos.dat
            self.write_dfpt_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
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
    
    # 创建mesh.conf  目的为了获得 total dos 
    def write_disp_mesh_conf(
        self,
        mesh_conf_dirpath, 
        species, 
        supercell,
        mp,
        ): 

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        mesh_conf_filepath = os.path.join(mesh_conf_dirpath, "mesh.conf")
        with open(mesh_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
    
    # 创建mesh.conf  目的为了获得 total dos  
    def write_dfpt_mesh_conf(
        self,
        mesh_conf_dirpath, 
        species, 
        supercell,
        mp,
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



