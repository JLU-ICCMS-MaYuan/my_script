import os
import re
import sys
import shutil
import logging
from pathlib import Path

from numpy import isclose

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
        pp_dir: str,
        **kwargs: dict,
        ):
        super(vasp_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path,
            pp_dir,
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
            self.kspacing = None
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

        if not hasattr(self, "lreal"):
            self.lreal = "Auto"
        
        if not hasattr(self, "symprec"):
            self.symprec = 1e-5

        if not hasattr(self, "mode"):
            self.mode = None 
        
        if not hasattr(self, "queue"):
            self.queue = None
            logger.debug("You didn't specify queue, so the program will not submit the job in any way")

        if not hasattr(self, "execmd"):
            logger.error("You must specify execute command, such as 'mpirun -np 48', 'bash', 'srun', 'srun --mpi=pmi2'")
            sys.exit(1)

        if not hasattr(self, "nbands"):
            self.nbands = None
            logger.warning("NBANDS had better to specify when you do eletronic calculation, espscially for COHP")
            logger.warning("Therefore, the default value will be determined by VASP automatically")

        # 关于并行计算的参数
        
        if not hasattr(self, "ncore"):
            self.ncore = None
            
        if not hasattr(self, "npar"):
            self.npar=None
        
        if not hasattr(self, "kpar"):
            self.kpar=None
            
        logger.debug(f"ncore = {self.ncore}  npar = {self.npar} kpar = {self.kpar}")

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

        # 关于电子态密度计算参数设置
        if not hasattr(self, "nedos"):
            self.nedos = 2000


    @classmethod
    def init_from_config1(cls, config: dict):

        work_path         = config['work_path']            ; del config['work_path']
        press             = config['press']                ; del config['press']
        submit_job_system = config['submit_job_system']    ; del config['submit_job_system']
        pp_dir            = config['pp_dir']               ; del config['pp_dir']
        input_file_path   = Path(config['input_file_path']); del config['input_file_path']
        self = cls(
            work_path=work_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            pp_dir=pp_dir,
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
        pp_dir: str,
        **kwargs: dict
        ):
        super(vasp_phonopara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            pp_dir,
            **kwargs,
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
            logger.debug("If you didn't set supercell, the default will be `supercell='1 1 1'` ")
            self.supercell = [1, 1, 1]

        if hasattr(self, "kdensity"):
            _kdensity = kwargs['kdensity'].split()
            self.kdensity = list(map(int, _kdensity))
            self.kspacing = None
        elif hasattr(self, "kspacing"):
            self.kdensity = None
            self.kspacing = kwargs["kspacing"]
        else:
            logger.debug("You didn't specify any information of k-mesh, so the program will set kdensity='40 40 40'")
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
                logger.debug("When mode=dfpt, the program will rename POSCAR to POSCAR-unit and copy SPOSCAR to POSCAR")
                shutil.move(poscar_file, poscar_init) 
                shutil.copy(sposcar_file, poscar_file)
            elif self.mode == 'disp':
                logger.debug("When mode=disp, the program will copy POSCAR to POSCAR-init and note that SPOSCAR won't be copied to POSCAR")
                shutil.copy(poscar_file, poscar_init)


class vasp_eletronpara(vasp_inputpara):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        pp_dir: str,
        **kwargs: dict,
        ):
        super(vasp_eletronpara, self).__init__(  # 这里特别重要，必须在super里面写vasp_eletronpara，而不是vasp_inputpara
            work_path, 
            press, 
            submit_job_system, 
            input_file_path,
            pp_dir,
            **kwargs, # 这里也特别重要，必须写**kwargs,因为vasp_eletronpara继承的函数vasp_inputpara有kwargs这个参数
            )
        
        self.set_default_inputpara(kwargs)
        try:
            self.mode = self.mode.split()
        except:
            logger.debug("You are using `vasp eletron module`, what you use mode is:")
            print(f"{self.mode}")
        
        if not hasattr(self, "vaspkitflag"):
            self.vaspkitflag = False
        else:
            self.vaspkitflag = eval(self.vaspkitflag)

        if not hasattr(self, "autoselect"):
            self.autoselect = False
        else:
            self.autoselect = eval(self.autoselect)


class vaspbatch_inputpara(vaspbatch_base, vasp_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: Path, 
        pp_dir: str,
        **kwargs: dict,
        ) -> None:

        super(vaspbatch_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            pp_dir, # 这里不需要加  **kwargs, 因为被继承的函数没有这个变量
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
        pp_dir: str,
        **kwargs: dict,
        ) -> None:

        super(vaspbatch_phonopara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path,
            pp_dir,
            )
        
        self.set_default_inputpara(kwargs)
        self.set_default_phonoinputpara(kwargs)


class vasp_mdpara(vasp_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        pp_dir: str,
        **kwargs: dict
        ):
        super(vasp_mdpara, self).__init__(
            work_path, 
            press,
            submit_job_system, 
            input_file_path, 
            pp_dir,
            **kwargs,
            )
        
        self.set_default_inputpara(kwargs)
        self.set_default_mdinputpara(kwargs)
    
    def set_default_mdinputpara(self):
         
        if not hasattr(self, "tebeg"):
            raise ValueError("You have to set tebeg(TEBEG)")

        if not hasattr(self, "teend"):
            raise ValueError("You have to set teend(TEEND))")

        if not hasattr(self, "potim"):
            self.potim = 0.5
        
        if not hasattr(self, "ibrion"):
            self.ibrion = 0
        elif hasattr(self, "ibrion") and int(self.ibrion) != 0:
            self.ibrion = 0

        elif hasattr(self, "potim") and isclose(float(self.potim), 0.1):
            raise("Your sure that POTIM = {} ??????".format(self.potim))
        else:
            self.potim = int(self.potim)

        if not hasattr(self, "nsw"):
            self.nsw = 20000
        else:
            self.nsw = int(self.nsw)

        if self.mode == 'nvt':
            self.isif = 2
            self.mdalgo = 2
            if not hasattr(self, "smass"):
                self.smass = 0 
        elif self.mode == 'npt':
            self.mdalgo = 3
            self.isif = 3
            if not hasattr(self, "langevin_gamma"):
                self.langevin_gamma = [10.0]*len(self.species)  
            else:
                self.langevin_gamma = list(map(float, self.langevin_gamma.split()))
            if not hasattr(self, "langevin_gamms_l"):
                self.langevin_gamms_l = 1 
            else:
                self.langevin_gamms_l = int(self.langevin_gamms_l)
            if not hasattr(self, "pmass"):
                self.pmass = 1.0
        elif self.mode == 'nve':
            self.mdalgo = 1
            if not hasattr(self, "smass"):
                self.smass = -3
            self.isif = 2
            self.andersen_prob = 0.0


if __name__ == "__main__":
    print(vaspbatch_inputpara.mro())