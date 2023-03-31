import os
import re
import logging
import shutil
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain


from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor


from vasp.config import config
from vasp.vasp_base import vasp_base
from vasp.vasp_inputpara import vasp_inputpara, vasp_phonopara, vaspbatch_inputpara, vaspbatch_phonopara
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
        _vasp_writeincar.opt_fine_incar(relax_inputpara.sub_workpath) 
    elif relax_inputpara.mode == 'rv3':
        _vasp_writeincar.opt_incar1(relax_inputpara.sub_workpath)
        _vasp_writeincar.opt_incar2(relax_inputpara.sub_workpath)
        _vasp_writeincar.opt_incar3(relax_inputpara.sub_workpath)
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
                _vasp_writeincar.opt_fine_incar(relax_inputpara.sub_workpath) 
            elif relax_inputpara.mode == 'rv3':
                _vasp_writeincar.opt_incar1(relax_inputpara.sub_workpath)
                _vasp_writeincar.opt_incar2(relax_inputpara.sub_workpath)
                _vasp_writeincar.opt_incar3(relax_inputpara.sub_workpath)
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
        _vasp_writeincar.disp_incar(phono_inputpara.sub_workpath)
    elif phono_inputpara.mode == 'dfpt':
        _vasp_writeincar.dfpt_incar(phono_inputpara.sub_workpath)

    # init the KPOINTS
    phono_inputpara.create_kpoints_by_pymatgen(
        phono_inputpara.sposcar_struct_type,
        phono_inputpara.sub_workpath.joinpath("KPOINTS"),
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
                _vasp_writeincar.disp_incar(phono_inputpara.sub_workpath)
            elif phono_inputpara.mode == 'dfpt':
                _vasp_writeincar.dfpt_incar(phono_inputpara.sub_workpath)
            # init the KPOINTS
            phono_inputpara.create_kpoints_by_pymatgen(
                phono_inputpara.sposcar_struct_type,
                phono_inputpara.sub_workpath.joinpath("KPOINTS"),
                phono_inputpara.kdensity,
                )
            # init the submit job script
            _vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(phono_inputpara)
            jobname = _vasp_writesubmit.write_submit_scripts()
            # submit the job
            _vasp_submitjob   = vasp_submitjob.init_from_phonoinput(phono_inputpara)
            if phono_inputpara.queue is not None:
                _vasp_submitjob.submit_mode2(jobname)


def vasp_properties(args):
        
    # read input para
    _config = config(args).read_config()

    # prepare the POSCAR POTCAR  
    properties_inputpara = vasp_inputpara.init_from_config1(_config)

    # init the INCAR
    _vasp_writeincar  = vasp_writeincar.init_from_properties(properties_inputpara)
    if properties_inputpara.mode == 'scf':
        _vasp_writeincar.scf_incar(
            properties_inputpara.sub_workpath
            )
        properties_inputpara.write_evenly_kpoints(
            properties_inputpara.kspacing, 
            properties_inputpara.sub_workpath
            )
        _vasp_writesubmit = vasp_writesubmit.init_from_properties(properties_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(properties_inputpara)
        if properties_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
    elif properties_inputpara.mode == 'eband':
        chgcar_src = properties_inputpara.work_path.joinpath("scf", "CHGCAR")
        chgcar_dst = properties_inputpara.work_path.joinpath("eband", "CHGCAR")
        shutil.copy(chgcar_src, chgcar_dst)
        properties_inputpara.write_highsymmetry_kpoints(
            properties_inputpara.ase_type, 
            properties_inputpara.sub_workpath.joinpath("KPOINTS"),
            )
        _vasp_writeincar.band_incar(properties_inputpara.sub_workpath)
        _vasp_writesubmit = vasp_writesubmit.init_from_properties(properties_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(properties_inputpara)
        if properties_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
        print("NOTES: ------------------------------ ")
        print("If you meet the erros in eband just like: ")
        print("    ERROR: charge density could not be read from file CHGCAR for ICHARG>10")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ in scf/OUTCAR are different from those in eband/OUTCAR ")
        print("    SOLUTION: You can let NGX,NGY,NGZ in eledos/INCAR be the same as scf/OUTCAR")
        print("If you meet the erros in eledos just like: ")
        print("    WARING: Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accuratecalculation. Thus, the results might be wrong. ")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ you'v customized aren't matched with the PREC=Accurate ")
        print("    SOLUTION: You can let PREC=Normal eband/INCAR")
    elif properties_inputpara.mode == 'eledos':
        chgcar_src = properties_inputpara.work_path.joinpath("scf", "CHGCAR")
        chgcar_dst = properties_inputpara.work_path.joinpath("eledos", "CHGCAR")
        shutil.copy(chgcar_src, chgcar_dst)
        print("NOTES:  ------------------------------ ")
        print("    KSPACING in `eledos` have to be twice than that in `scf` ")
        properties_inputpara.write_evenly_kpoints(
            properties_inputpara.kspacing, 
            properties_inputpara.sub_workpath
            )
        _vasp_writeincar.eledos_incar(properties_inputpara.sub_workpath)
        _vasp_writesubmit = vasp_writesubmit.init_from_properties(properties_inputpara)
        jobname = _vasp_writesubmit.write_submit_scripts()
        _vasp_submitjob = vasp_submitjob.init_from_relaxinput(properties_inputpara)
        if properties_inputpara.queue is not None:
            _vasp_submitjob.submit_mode1(jobname)
        print("NOTES: ------------------------------ ")
        print("if you meet the erros just like: ")
        print("    ERROR: charge density could not be read from file CHGCAR for ICHARG>10")
        print("    Possible reason is that NGX, NGY, NGZ in scf/OUTCAR are different from those in eledos/OUTCAR ")
        print("    If you wanna to solve it, you need let NGX,NGY,NGZ in eledos/INCAR be the same as scf/OUTCAR")
        print("If you meet the erros just like: ")
        print("    WARING: Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accuratecalculation. Thus, the results might be wrong. ")
        print("    ANALYSIS: Possible reason is that NGX, NGY, NGZ you'v customized aren't matched with the PREC=Accurate ")
        print("    SOLUTION: You can let PREC=Normal in eledos/INCAR ")
    else:
        print(vasp_writeincar.mode)


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
    
    # 绘制 phonoband 
    def post_progress_phono_band(self):

        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            print("WARNING: ----------------------------- ")
            raise ValueError("you have to specify the supercell=[?,?,?]. If you didn't specify it, maybe somthing wrong will occur !")

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [8,8,8]
            print("NOTES: ------------------------------ ")
            print("you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if self.mode == "dispband":

            _disp_num = len(list(Path(self.work_path).glob("disp-*")))
            disp_num = str(_disp_num).rjust(3, '0')

            cwd = os.getcwd()
            os.chdir(self.work_path)
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
            os.system("phonopy -p -s band.conf -c POSCAR-init")
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)

        elif self.mode == "dfptband":

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
            os.system("phonopy --dim='{}' -p -s band.conf -c POSCAR-init".format(' '.join(list(map(str, self.supercell)))))
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)
    
    # 绘制 phonodos 
    def post_progress_phono_dos(self):

        if "phdos" in self._config:
            self.pdos = self._config['pdos']
        else:
            self.pdos = "AUTO"

        if self.mode == "dispprog": 
            # 获得total_dos.dat
            self.write_disp_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf") # -p: dos plot   -t: thermal properties print
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

        elif self.mode == "dfptprog":
            # 获得total_dos.dat
            self.write_dfpt_mesh_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                self.mp,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf") # -p: dos plot   -t: thermal properties print
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
        from pymatgen.electronic_structure.plotter import BSDOSPlotter

        vasprunxml_path = self.work_path.joinpath("vasprun.xml")
        vasprun = Vasprun(vasprunxml_path, parse_projected_eigen=True)

        eband = vasprun.get_band_structure(line_mode=True)

        e_fermi_fromvasp = vasprun.efermi
        e_fermi_frompmg  = vasprun.calculate_efermi(float=0.0001)

        # set figure parameters, draw figure
        eband_fig = BSDOSPlotter(bs_projection=None, dos_projection=None, vb_energy_range=5, fixed_cb_energy=5)
        eband_fig.get_plot(bs=eband)
        ebandpng_path = self.work_path.joinpath('eband.png')
        plt.savefig(ebandpng_path, img_format='png')
    
    # 绘制 eledos
    def post_progress_eletron_dos(self):

        import matplotlib.pyplot as plt
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import BSDOSPlotter

        vasprunxml_path = self.work_path.joinpath("vasprun.xml")
        vasprun = Vasprun(vasprunxml_path)

        eledos = vasprun.complete_dos

        e_fermi_fromvasp = vasprun.efermi
        e_fermi_frompmg  = vasprun.calculate_efermi(tol=0.0001)
        print("NOTES: ------------------------------ ")
        print("    You have to confirm that the Fermi energy is from scf/OUTCAR. Because the Fermi energy in dos/DOSCAR is not accurate")
        print("    You cat use 'grep E-fermi scf/OUTCAR' to check the Fermi energy! ")
        # set figure parameters, draw figure
        eledos_fig = BSDOSPlotter(bs_projection=None, dos_projection=None, vb_energy_range=5, fixed_cb_energy=5)
        eledos_fig.get_plot(bs=eledos)
        eledospng_path = self.work_path.joinpath('eband.png')
        plt.savefig(eledospng_path, img_format='png')
    
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



