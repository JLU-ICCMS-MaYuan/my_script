import os
import re
import sys
import shutil
import logging
from pathlib import Path

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

from vasp.vasp_base import vasp_base, vaspbatch_base

logger = logging.getLogger(__name__)

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
        mode: str,
        **kwargs: dict,
        ):
        super(vasp_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path,
            mode,
            )

        self.set_default_inputpara(kwargs)

    def set_default_inputpara(self, kwargs):

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
        else:
            self.kspacing = float(self.kspacing)

        if not hasattr(self, "ismear"):
            self.ismear = 0

        if not hasattr(self, "sigma"):
            self.sigma = 0.01

        if not hasattr(self, "ediff"):
            self.ediff = 1E-07

        if not hasattr(self, "ediffg"):
            self.ediffg = -0.01

        if not hasattr(self, "ibrion"):
            self.ibrion = 2

        if not hasattr(self, "isif"):
            self.isif = 3

        if not hasattr(self, "potim"):
            self.potim = 0.1

        if not hasattr(self, "nelm"):
            self.nelm = 200

        if not hasattr(self, "ncore"):
            self.ncore = 4

        if not hasattr(self, "lreal"):
            self.lreal = "Auto"
        
        if not hasattr(self, "symprec"):
            self.symprec = 1e-5

        if not hasattr(self, "mode"):
            self.mode = None 
        
        if not hasattr(self, "queue"):
            self.queue = None
            print("    You didn't specify queue, so the program will not submit the job in any way")

        if not hasattr(self, "core"):
            print("Error: ----------------------")
            print("    You must specify the number of core, such as 'core=48'")
            print("-----------------------------")
            sys.exit(1)

        if not hasattr(self, "nbands"):
            self.nbands = 100
            print("Warning: --------------------")
            print("    NBANDS had better to specify when you do eletronic properties calculation!!!")
            print("    Therefore, the default value,  NBANDS = 100")
            print("-----------------------------")

        # 关于并行计算的参数
        if not hasattr(self, "npar"):
            self.npar=4
        
        # 关于磁性的参数设置：
        if not hasattr(self, "isym"):
            self.isym = 2
        
        if not hasattr(self, "ispin"):
            self.ispin = 1
        
        if not hasattr(self, "magmom"):
            self.magmom = None

        if not hasattr(self, "lorbit"):
            self.lorbit = 11
        
        if not hasattr(self, "lasph"):
            self.lasph = ".TRUE."
        
        if not hasattr(self, "gga"):
            self.gga = "PS"

        # 关于DFT+U+J的参数设置
        if not hasattr(self, "ldau"):
            self.ldau = None

        if not hasattr(self, "ldautype"):
            self.ldautype = 2

        if not hasattr(self, "ldaul"):
            self.ldaul = None

        if not hasattr(self, "ldauu"):
            self.ldauu = None

        if not hasattr(self, "ldauj"):
            self.ldauj = None

        if not hasattr(self, "lmaxmix"):
            self.lmaxmix = 6

    @classmethod
    def init_from_config1(cls, config: dict):

        work_path         = config['work_path']            ; del config['work_path']
        press             = config['press']                ; del config['press']
        submit_job_system = config['submit_job_system']    ; del config['submit_job_system']
        input_file_path   = Path(config['input_file_path']); del config['input_file_path']
        mode              = config['mode']                 ; del config['mode']
        self = cls(
            work_path=work_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            mode=mode,
            **config,
        )
        return self


class vasp_phonopara(vasp_inputpara):
    """
    add more paramaters for vasp phono calculate
    It inherits the vasp_inputpara class.
    
    \nNote:
        please attention to specify the lreal(LREAL) and ncore(NCORE) 可能会影响计算速度
        The formate is:     lreal=.TRUE. or lreal=.FALSE.
                            ncore=4      or ncore=1 
    """

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        mode: str, 
        **kwargs: dict
        ):
        super(vasp_phonopara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            mode, 
            **kwargs
            )
        
        self.set_default_inputpara(kwargs)
        self.set_default_phonoinputpara(kwargs)

    def set_default_phonoinputpara(self, kwargs):

        if hasattr(self, "kpoints"):
            _kpoints = kwargs['kpoints'].split()
            self.kpoints = list(map(int, _kpoints))
        else:
            self.kpoints = [None, None, None]

        if hasattr(self, "supercell"):
            _supercell = re.findall(r"\d+", kwargs['supercell'])
            self.supercell = list(map(int, _supercell))
        else:
            print("If you didn't set supercell, the default will be `supercell='1 1 1'` ")
            self.supercell = [1, 1, 1]

        if hasattr(self, "kdensity"):
            _kdensity = kwargs['kdensity'].split()
            self.kdensity = list(map(int, _kdensity))
            self.kspacing = None
        elif hasattr(self, "kspacing"):
            self.kdensity = None
            self.kspacing = kwargs["kspacing"]
        else:
            print("You didn't specify any information of k-mesh, so the program will set kdensity='40 40 40'")
            self.kdensity = [40, 40, 40]
            self.kspacing = None
        
        if self.mode == "disp" or self.mode == "dfpt":
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system('phonopy -d --dim="{} {} {}"'.format(
                        self.supercell[0],
                        self.supercell[1],
                        self.supercell[2],
                    ))
            os.chdir(cwd)

            poscar_file = Path(self.work_path).joinpath("POSCAR")
            poscar_init = Path(self.work_path).joinpath("POSCAR-init")
            sposcar_file= Path(self.work_path).joinpath("SPOSCAR")
            self.sposcar_ase_type    = read(sposcar_file)
            self.sposcar_struct_type = AseAtomsAdaptor.get_structure(self.sposcar_ase_type) 
            if self.mode == 'dfpt':
                print("When mode=dfpt, the program will rename POSCAR to POSCAR-unit and copy SPOSCAR to POSCAR")
                shutil.move(poscar_file, poscar_init) 
                shutil.copy(sposcar_file, poscar_file)
            elif self.mode == 'disp':
                print("When mode=disp, the program will copy POSCAR to POSCAR-init and note that SPOSCAR won't be copied to POSCAR")
                shutil.copy(poscar_file, poscar_init)


class vasp_eletronpara(vasp_inputpara):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        mode: str,
        **kwargs: dict,
        ):
        super(vasp_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path,
            mode,
            )
        
        self.set_default_inputpara(kwargs)
        try:
            self.mode = self.mode.split()
        except:
            print("You are using `vasp eletron module`, what you use mode is:")
            print(f"{self.mode}")
            

class vaspbatch_inputpara(vaspbatch_base, vasp_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: Path, 
        mode: str,
        **kwargs: dict,
        ) -> None:

        super(vaspbatch_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, mode
            )
        
        '''
            由于第一次使用多继承, 这里做一点简短的说明
            根据多继承的mro链顺序, 得到vaspbatch_inputpara类的mro顺序为 : 
            [
                <class '__main__.vaspbatch_inputpara'>, 
                <class 'vasp.vasp_base.vaspbatch_base'>, 
                <class '__main__.vasp_inputpara'>, 
                <class 'vasp.vasp_base.vasp_base'>,
                <class 'object'>
            ]
            所有在使用super调用继承init方法时, 会优先寻找vaspbatch_inputpara中的init方法, 
            如果vaspbatch_inputpara中的init方法不存在, 则会寻找vasp_inputpara的init方法.
            但是vaspbatch_inputpara类继承的是vasp_base类。那么vaspbatch_inputpara类会不会继承vaspbatch_inputpara类的父类vasp_base类呢？
            答应是 : 不会。
            因为python3的新类继承遵循的原则是 : 广度优先。不会追根溯源的继承最根本的父类的。所有vaspbatch_inputpara继承到vaspbatch_inputpara就停止了.
        '''

        self.set_default_inputpara(kwargs)


class vaspbatch_phonopara(vaspbatch_base, vasp_phonopara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: Path, 
        mode: str,
        **kwargs: dict,
        ) -> None:

        super(vaspbatch_phonopara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, mode
            )
        
        self.set_default_inputpara(kwargs)
        self.set_default_phonoinputpara(kwargs)


if __name__ == "__main__":
    print(vaspbatch_inputpara.mro())