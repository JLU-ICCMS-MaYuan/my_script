import os
import logging
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain

from config import config
from vasp_inputpara import vasp_inputpara 
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
        self.phono_inputpara  = vasp_inputpara.init_from_config1(self._config)

        # init the INCAR
        self.vasp_writeincar  = vasp_writeincar.init_from_phonoinput(self.phono_inputpara)

        # init the KPOINTS
        if self.phono_inputpara.kpoints != [None, None, None]:
            self.vasp_kpoints     = vasp_writekpoints.init_from_inputpara(self.phono_inputpara)

        # init the submit job script
        self.vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(self.phono_inputpara)

        # submit the job
        self.vasp_submitjob   = vasp_submitjob.init_from_phonoinput(self.phono_inputpara)

        if  self.phono_inputpara.mode == "dispprog" or\
            self.phono_inputpara.mode == "dfptprog":
            self.post_progress(
                self.phono_inputpara.mode
            )

    def post_progress(
        self, 
        mode,
        ):

        if mode == "dispprog":

            _disp_num = len(list(Path(self.phono_inputpara.work_underpressure).glob("disp-*")))
            disp_num = str(_disp_num).rjust(3, '0')

            cwd = os.getcwd()
            os.chdir(self.phono_inputpara.work_underpressure)
            os.system("phonopy -f disp-{001..%s}/vasprun.xml" %(disp_num))
            os.chdir(cwd) 

            self.write_Disp_band_conf(
                self.phono_inputpara.work_underpressure, 
                self.phono_inputpara.species, 
                self.phono_inputpara.supercell, 
                special_points,
                path_coords
                )     
        elif mode == "dfptprog":

            cwd = os.getcwd()
            os.chdir(self.phono_inputpara.work_underpressure)
            os.system("phonopy --fc vasprun.xml")
            os.chdir(cwd)
            special_points, path_coords = self.get_hspp()
            
            self.write_dfpt_band_conf(
                self.phono_inputpara.work_underpressure, 
                self.phono_inputpara.species, 
                self.phono_inputpara.supercell, 
                special_points,
                path_coords
                )     

            cwd = os.getcwd()
            os.chdir(self.phono_inputpara.work_underpressure)
            os.system("phonopy -p -s band.conf -c POSCAR-init")
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)
    
    def get_hspp(self):

        ltype   = self.phono_inputpara.ase_type.cell.get_bravais_lattice()
        pstring = ltype.special_path
        _plist  = [[ p for p in pp] for pp in pstring.split(",")]

        logger.info(f"the high symmetry points path is {_plist}")

        high_symmetry_type = input(
            "please input the mode you want\n",
            "'all_points'\n",
            "'main_points'\n",
            "Nothing to input\n")
        if not high_symmetry_type:
            high_symmetry_type = "all_points" # default
        if "," in pstring:
            if high_symmetry_type == "all_points":
                plist = list(chain.from_iterable(_plist))
                logger.info(f"the choosed high symmetry points path is \n {plist}")
            elif high_symmetry_type == "main_points":
                plist = max(_plist, key=len)
                logger.info(f"the choosed high symmetry points path is \n {plist}")
        else:
            plist = [ pp for pp in pstring]
            logger.info(f"the choosed high symmetry points path is \n {plist}")

        special_points   = ltype.get_special_points()
        path_coords      = [list(special_points[pname]) for pname in plist]

        return special_points, path_coords
        
    def write_disp_band_conf(
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
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("NPOINTS=101             \n")
            f.write("EIGENVECTORS=.TRUE.     \n")
            f.write("BAND_LABELS={}          \n".format(' '.join(special_points)))
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
                self.relax_inputpara  = vasp_inputpara(
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
                if self.phono_inputpara.kpoints != [None, None, None]:
                    self.vasp_kpoints     = vasp_writekpoints.init_from_inputpara(self.phono_inputpara)

                # init the submit job script
                self.vasp_writesubmit = vasp_writesubmit.init_from_phonoinput(self.phono_inputpara)

                # submit the job
                self.vasp_submitjob   = vasp_submitjob.init_from_phonoinput(self.phono_inputpara)

                if  self.phono_inputpara.mode == "dispprog" or\
                    self.phono_inputpara.mode == "dfptprog":
                    self.post_progress(
                        self.phono_inputpara.mode
                    )