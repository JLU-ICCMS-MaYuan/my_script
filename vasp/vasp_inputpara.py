import os
import re
import shutil
from pathlib import Path
from argparse import ArgumentParser

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

from vasp_base import vaspbase

class vasp_inputpara(vaspbase):

    def __init__(
        self, 
        args: ArgumentParser,
        ) -> None:
        super().__init__(args)

        if 'encut' in self._config:
            self.encut = self._config['encut']
        else:
            self.encut = 600
        
        if 'kspacing' in self._config:
            self.kspacing = self._config['kspacing']
        else:
            self.kspacing = 0.3

        if 'kpoints' in self._config:
            _kpoints = re.findall(r"\d+", self._config['kpoints'])
            self.kpoints = list(map(int, _kpoints))
        else:
            self.kpoints = [None, None, None]

        if 'ismear' in self._config:
            self.ismear = self._config['ismear']
        else:
            self.ismear = 0

        if 'sigma' in self._config:
            self.sigma = self._config['sigma']
        else:
            self.sigma = 0.02

        if 'ediff' in self._config:
            self.ediff = self._config['ediff']
        else:
            self.ediff = 1e-8

        if 'ibrion' in self._config:
            self.ibrion = self._config['ibrion']
        else:
            self.ibrion = 2

        if 'isif' in self._config:
            self.isif = self._config['isif']
        else:
            self.isif = 3

        if 'optim' in self._config:
            self.optim = self._config['optim']
        else:
            self.optim = 0.1
        
        if 'mode'  in self._config:
            self.mode = self._config['mode']
        else:
            self.mode = None 
        
        if 'supercell' in self._config:
            _supercell = re.findall(r"\d+", self._config['supercell'])
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
            shutil.copy(poscar_file, poscar_init) 
        
        if self.mode == "dispprog":
            pass

        if self.mode == "dfptprog":
            pass
