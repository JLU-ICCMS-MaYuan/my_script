import os
import re
import shutil
from pathlib import Path
from argparse import ArgumentParser
from turtle import pos

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

from config import config
from vasp_base import vasp_base 

class vasp_inputpara(vasp_base):
    """
    Preparing the Input File: POSCAR, POTCAR and many parameters of INCAR for single structure.
    It inherits the vaspBase class.
    It can't prepare the input file for batch computing.
    """
    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict,
        ):
        super(vasp_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path
            )

        for key, value in kwargs.items():
            if key != "work_path" and \
               key != "press" and \
               key != "submit_job_system" and \
               key !="input_file_path":
                setattr(self, key, value)

        if not hasattr(self, "encut"):
            self.encut = 600
        
        if not hasattr(self, "kspacing"):
            self.kspacing = 0.3

        if hasattr(self, "kpoints"):
            _kpoints = re.findall(r"\d+", kwargs['kpoints'])
            self.kpoints = list(map(int, _kpoints))
        else:
            self.kpoints = [None, None, None]

        if not hasattr(self, "ismear"):
            self.ismear = 0

        if not hasattr(self, "sigma"):
            self.sigma = 0.02

        if not hasattr(self, "ediff"):
            self.ediff = 1e-8

        if not hasattr(self, "ibrion"):
            self.ibrion = 2

        if not hasattr(self, "isif"):
            self.isif = 3

        if not hasattr(self, "potim"):
            self.potim = 0.1
        
        if not hasattr(self, "mode"):
            self.mode = None 
        
        if not hasattr(self, "queue"):
            self.queue = "xieyu"

        if hasattr(self, "supercell"):
            _supercell = re.findall(r"\d+", kwargs['supercell'])
            self.supercell = list(map(int, _supercell))
        else:
            self.supercell = [1, 1, 1]
        
        if self.mode == "disp" or self.mode == "dfpt":
            cwd = os.getcwd()
            os.chdir(self.work_underpressure)
            os.system('phonopy -d --dim="{} {} {}"'.format(
                        self.supercell[0],
                        self.supercell[1],
                        self.supercell[2],
                    ))
            os.chdir(cwd)

            poscar_file = Path(self.work_underpressure).joinpath("POSCAR")
            poscar_init = Path(self.work_underpressure).joinpath("POSCAR-init")
            sposcar_file= Path(self.work_underpressure).joinpath("SPOSCAR")
            self.sposcar_ase_type    = read(sposcar_file)
            self.sposcar_struct_type = AseAtomsAdaptor.get_structure(self.sposcar_ase_type) 
            if self.mode == 'dfpt':
                shutil.move(poscar_file, poscar_init) 
                shutil.copy(sposcar_file, poscar_file)
            elif self.mode == 'disp':
                shutil.copy(poscar_file, poscar_init)
                

        if self.mode == "dispprog":
            pass

        if self.mode == "dfptprog":
            pass

    @classmethod
    def init_from_config1(cls, config: dict):

        work_path         = config['work_path']            ; del config['work_path']
        press             = config['press']                ; del config['press']
        submit_job_system = config['submit_job_system']    ; del config['submit_job_system']
        input_file_path   = Path(config['input_file_path']); del config['input_file_path']
        self = cls(
            work_path=work_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            **config,
        )
        return self

    @classmethod
    def init_from_config2(cls, input_file_path, config: dict):

        self = cls(
            work_path=input_file_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            **config,
        )
        return self