import os
import re
import logging
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

from config import config
from vasp_base import vasp_base
from vasp_inputpara import vasp_inputpara, vasp_phonopara 
from vasp_writeincar import vasp_writeincar
from vasp_writekpoints import vasp_writekpoints
from vasp_writesubmit import vasp_writesubmit
from vasp_submitjob import vasp_submitjob

logger = logging.getLogger("vasp_run")

class vasp_relax:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()

        # prepare the POSCAR POTCAR  
        self.relax_inputpara  = vasp_inputpara.init_from_config1(self._config)

        # init the INCAR
        self.vasp_writeincar  = vasp_writeincar.init_from_relaxinput(self.relax_inputpara)

        # init the submit job script
        self.vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(self.relax_inputpara)

        # submit the job
        self.vasp_submitjob   = vasp_submitjob.init_from_relaxinput(self.relax_inputpara)


class vaspbatch_relax:

    def __init__(self, args: ArgumentParser) -> None:

        self._config = config(args).read_config()
        self.input_dir_path = Path(self._config['input_file_path'])
        if self.input_dir_path.is_dir():
            self.input_files_path = list(self.input_dir_path.glob("*.vasp"))
            work_path         = self._config['work_path']        ; del self._config['work_path']
            press             = self._config['press']            ; del self._config['press']
            submit_job_system = self._config['submit_job_system']; del self._config['submit_job_system']
            pass                                                 ; del self._config['input_file_path']
            mode              = self._config['mode']             ; del self._config['mode']
            for input_file_path in self.input_files_path:
                # prepare the POSCAR POTCAR  
                self.relax_inputpara  = vasp_inputpara(
                    work_path=work_path,
                    press=press,
                    submit_job_system=submit_job_system,
                    input_file_path=input_file_path,
                    mode=mode,
                    **self._config
                    )

                # init the INCAR
                self.vasp_writeincar  = vasp_writeincar.init_from_relaxinput(self.relax_inputpara)

                # init the submit job script
                self.vasp_writesubmit = vasp_writesubmit.init_from_relaxinput(self.relax_inputpara)

                # submit the job
                self.vasp_submitjob   = vasp_submitjob.init_from_relaxinput(self.relax_inputpara)


class vasp_phono:

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()
        
        # prepare the POSCAR POTCAR  
        self.phono_inputpara  = vasp_phonopara.init_from_config1(self._config)

        # init the INCAR
        self.vasp_writeincar  = vasp_writeincar.init_from_phonoinput(self.phono_inputpara)

        # init the KPOINTS
        self.vasp_kpoints     = vasp_writekpoints.init_from_inputpara(self.phono_inputpara)

        # init the submit job script
        self.vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(self.phono_inputpara)

        # submit the job
        self.vasp_submitjob   = vasp_submitjob.init_from_phonoinput(self.phono_inputpara)


class vaspbatch_phono(vasp_phono):

    def __init__(self, args: ArgumentParser) -> None:

        self._config = config(args).read_config()
        self.input_dir_path = Path(self._config['input_file_path'])
        if self.input_dir_path.is_dir():
            self.input_files_path = list(self.input_dir_path.glob("*.vasp"))
            work_path         = self._config['work_path']        ; del self._config['work_path']
            press             = self._config['press']            ; del self._config['press']
            submit_job_system = self._config['submit_job_system']; del self._config['submit_job_system']
            pass                                                 ; del self._config['input_file_path']
            mode              = self._config['mode']             ; del self._config['mode']
            for input_file_path in self.input_files_path:
                # prepare the POSCAR POTCAR  
                self.phono_inputpara  = vasp_phonopara(
                    work_path=work_path,
                    press=press,
                    submit_job_system=submit_job_system,
                    input_file_path=input_file_path,
                    mode=mode,
                    **self._config
                    )

                # init the INCAR
                self.vasp_writeincar  = vasp_writeincar.init_from_phonoinput(self.phono_inputpara)

                # init the KPOINTS
                self.vasp_kpoints     = vasp_writekpoints.init_from_inputpara(self.phono_inputpara)

                # init the submit job script
                self.vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(self.phono_inputpara)

                # submit the job
                self.vasp_submitjob   = vasp_submitjob.init_from_phonoinput(self.phono_inputpara)


class vasp_processdata(vasp_base):

    def __init__(self, args: ArgumentParser) -> None:

        # read input para
        self._config = config(args).read_config()
        self.input_file_path = Path(self._config['input_file_path'])
        self.work_path       = Path(self._config['work_path'])
        self.mode            = self._config['mode']
        # the default value of spectrum and dos is False, 
        # that means the program will not create the phono-spectrum, total-dos, pdos
        self.spectrum        = eval(self._config.get('spectrum', "False"))
        self.dos             = eval(self._config.get('dos', "False"))

        self.ase_type          = read(self.input_file_path)
        self.struct_type       = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)
        
        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            raise ValueError("you have to specify the supercell=[?,?,?]")

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [8,8,8]
            logger.warning("you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if "pdos" in self._config:
            self.pdos = self._config['pdos']
        else:
            self.pdos = "AUTO"

        if self.spectrum:
            self.post_progress_phono_spectrum()
        if self.dos:
            self.post_progress_phono_dos()



    def post_progress_phono_spectrum(self):

        if self.mode == "dispprog":

            _disp_num = len(list(Path(self.work_path).glob("disp-*")))
            disp_num = str(_disp_num).rjust(3, '0')

            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -f disp-{001..%s}/vasprun.xml" %(disp_num))
            os.chdir(cwd) 
            path_name_list, path_coords = self.get_hspp()

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

        elif self.mode == "dfptprog":

            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy --fc vasprun.xml")
            os.chdir(cwd)
            special_points, path_coords = self.get_hspp()
            
            self.write_dfpt_band_conf(
                self.work_path, 
                self.species, 
                self.supercell, 
                special_points,
                path_coords
                )     

            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -s band.conf -c POSCAR-init")
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)

    def get_hspp(self):

        ltype   = self.ase_type.cell.get_bravais_lattice()
        pstring = ltype.special_path
        _plist  = [[ p for p in pp] for pp in pstring.split(",")]

        logger.info(f"the high symmetry points path is {_plist}")

        print(
            "please input the mode you want\n",
            "'all_points'\n",
            "'main_points'\n",
            "Nothing to input\n"
            )
        high_symmetry_type = input()

        if not high_symmetry_type:
            high_symmetry_type = "all_points" # default
        if "," in pstring:
            if high_symmetry_type == "all_points":
                path_name_list = list(chain.from_iterable(_plist))
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
            elif high_symmetry_type == "main_points":
                path_name_list = max(_plist, key=len)
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in pstring]
            logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
        
        # 这里的高对称点是指：这个倒空间的布里渊区有几种高对称点，并不是某一个路径的高对称点的列表。
        # 如果想获得某一个路径下高对称点的坐标，需要按照path_name_list的顺序依次获得相应的坐标。
        special_points   = ltype.get_special_points()
        path_coords      = [list(special_points[pname]) for pname in path_name_list]

        return path_name_list, path_coords
        
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



    def post_progress_phono_dos(self):

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
            self.write_disp_pdos_conf(
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
            self.write_dfpt_pdos_conf(
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
    def write_disp_pdos_conf(
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
    def write_dfpt_pdos_conf(
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