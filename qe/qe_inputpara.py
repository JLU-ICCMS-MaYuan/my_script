import os
import re
import sys
import shutil
import logging
from pathlib import Path
from itertools import chain
from math import ceil

import numpy as np

from qe.qe_base import qe_base

logging.basicConfig(
    level = logging.INFO, 
    format='%(asctime)s | %(name)s | %(levelname)s ---- %(message)s'
    )
logger = logging.getLogger(__name__)

class qe_inputpara(qe_base):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict,
        ):
        super(qe_inputpara, self).__init__(
            work_path,
            press,
            submit_job_system,
            input_file_path
        )

        for key, value in kwargs.items():
            setattr(self, key, value)

        if not hasattr(self, "qe_workflow"):
            raise AttributeError("there is no attribution of qe_workflow")

        if not hasattr(self, "mode"):
            raise AttributeError("there is no attribution of mode")

        if not hasattr(self, "core"):
            raise AttributeError("there is no attribution of core")

        if not hasattr(self, "npool"):
            self.npool = 1
            print("The program will set npool=1 in your submitjob scripts")
 
        if not hasattr(self, "queue"):
            self.queue = None
            print("You didn't specify queue, so the program will not submit the job in any way")

        print('\n'); print("if you want to run `relax` `scffit` `fit`, you had better set these values!")
        # &CONTROL
        if not hasattr(self, "forc_conv_thr"):
            self.forc_conv_thr = "1.0d-5"
            print("You didn't set the `forc_conv_thr` ! The program will use default value: forc_conv_thr=1.0d-5")

        if not hasattr(self, "etot_conv_thr"):
            self.etot_conv_thr = "1.0d-7"
            print("You didn't set the `etot_conv_thr` ! The program will use default value: etot_conv_thr=1.0d-7")

        # &SYSTEM
        if not hasattr(self, "occupations"):
            self.occupations = "smearing"
            print("You didn't set the `occupations` !   The program will use default value: occupations=smearing")

        if not hasattr(self, "smearing"):
            self.smearing = "gauss"
            print("You didn't set the `smearing` !      The program will use default value: smearing=gauss")

        # self.smearing = methfessel-paxton 做scffit 和 scf 时候用这个参数
        if not hasattr(self, "degauss"):
            self.degauss = "0.005"
            print("You didn't set the `degauss` !       The program will use default value: degauss=0.005")

        if not hasattr(self, "ecutwfc"):
            self.ecutwfc = "60"
            print("You didn't set the `ecutwfc` !       The program will use default value: ecutwfc=60")

        if not hasattr(self, "ecutrho"):
            self.ecutrho = "720"
            print("You didn't set the `ecutrho` !       The program will use default value: ecutrho=720")

        if not hasattr(self, "lspinorb"):
            self.lspinorb = "false"
            print("You didn't set the `lspinorb` !      The program will use default value: lspinorb=false")
        else:
            print("Please carefully check the bool value of `lspinorb` you just set. Its format must be `false` or `true` without capital")

        if not hasattr(self, "noncolin"):
            self.noncolin = "false"
            print("You didn't set the `noncolin` !      The program will use default value: noncolin=false")
        elif self.lspinorb == "true":
            self.noncolin = "true"
            print("Because lspinorb = true, so the noncolin=true")
        else:
            print("Please carefully check the bool value of `noncolin` you just set. Its format must be `false` or `true` without capital")

        if not hasattr(self, "la2F"):
            self.la2F = "true"
            print("You didn't set the `la2F` !          The program will use default value: la2F=true. ")
            print("But in relax mode and scf mode, it doesn't exist ! It only exist in scffit mode")
        else:
            print("Please carefully check the bool value of `la2F` you just set. Its format must be `false` or `true` without capital")

        # &ELECTRONS
        if not hasattr(self, "diagonalization"):
            self.diagonalization = "david"
            print("You didn't set the `diagonalization`! The program will use default value: diagonalization=david")
        
        if not hasattr(self, "conv_thr"):
            self.conv_thr = "1.0d-8"
            print("You didn't set the `conv_thr` !      The program will use default value: conv_thr=1.0d-8")
            #  做结构弛豫 1.0-d8
            #  做scffit 和 scf 时候用1.0d-9
        if not hasattr(self, "mixing_beta"):
            self.mixing_beta = "0.7"
            print("You didn't set the `mixing_beta` !   The program will use default value: mixing_beta=0.7")
            #  做scffit 和 scf 时候用0.8
        if not hasattr(self, "electron_maxstep"):
            self.electron_maxstep = "200"
            print("You didn't set the `electron_maxstep`! The program will use default value: electron_maxstep=200")

        # &CELL
        if not hasattr(self, "press_conv_thr"):
            self.press_conv_thr = "0.01"
            print("You didn't set the `press_conv_thr`! The program will use default value: press_conv_thr=0.01")

        # &kpoints
        if hasattr(self, "kpoints_dense"):
            _kpoints_dense = self.kpoints_dense.split()
            self.kpoints_dense = list(map(int, _kpoints_dense))
        else:
            self.kpoints_dense = [16, 16, 16]    

        if hasattr(self, "kpoints_sparse"):
            _kpoints_sparse = self.kpoints_sparse.split()
            self.kpoints_sparse = list(map(int, _kpoints_sparse))
        else:
            self.kpoints_sparse = [8, 8, 8]    

    @classmethod
    def init_from_config(cls, config: dict):

        work_path         = config['work_path']            ; del config['work_path']
        press             = config['press']                ; del config['press']
        submit_job_system = config['submit_job_system']    ; del config['submit_job_system']
        input_file_path   = config['input_file_path']      ; del config['input_file_path']
        
        self = cls(
            work_path=work_path,
            press=press,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            **config,
        )
        return self


class qephono_inputpara(qe_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict
        ):
        super(qephono_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs
            )
        
        print('\n'); print("if you want to run `phono`, you had better set these values!")
        dyn0_names = Path(self.work_underpressure).joinpath(f"{self.system_name}.dyn0")
        if hasattr(self, "qpoints"):
            _qpoints = self.qpoints.split()
            self.qpoints = list(map(int, _qpoints))
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]
        elif dyn0_names.exists():
            print(f"You didn't set the `qpoints` !        The program will qpoints in read {self.system_name}.dyn0 file")
            self.qpoints = self.get_qpoints(dyn0_path=dyn0_names)
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints] 
        else:
            self.qpoints = [4,4,4]
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]

        if not hasattr(self, "tr2_ph"):
            self.tr2_ph = "1.0d-16"
            print(f"You didn't set the `tr2_ph` !          The program will use default value: tr2_ph=1.0d-16")

        if not hasattr(self, "electron_phonon"):
            self.electron_phonon="interpolated"
            print(f"You didn't set the `electron_phonon` ! The program will use default value: electron_phonon=interpolated")

        if not hasattr(self, "el_ph_nsigma"):
            self.el_ph_nsigma = "10"
            print(f"You didn't set the `el_ph_nsigma` !    The program will use default value: el_ph_nsigma=10")
            
        if not hasattr(self, "el_ph_sigma"):
            self.el_ph_sigma = "0.005"
            print(f"You didn't set the `el_ph_sigma` !     The program will use default value: el_ph_sigma=0.005")
        
        if not hasattr(self, "alpha_mix"):
            self.alpha_mix = "0.5"
            print(f"You didn't set the `alpha_mix` !       The program will use default value: alpha_mix=0.5")
        
        if not hasattr(self, "dyn0_flag"):
            self.dyn0_flag = False
        
        else:
            self.dyn0_flag = eval(self.dyn0_flag)

        dyn0_names = Path(self.work_underpressure).joinpath(f"{self.system_name}.dyn0")
        if dyn0_names.exists():
            self.qtot, self.qirreduced, self.qirreduced_coords= self.get_q_from_dyn0(dyn0_path=dyn0_names) 
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            print(f"{self.system_name}.dyn0 doesn't exist in {self.work_underpressure}. The qtot, qirreduced, qirreduced_coords will not get values.")

        if not hasattr(self, "qtot"):
            self.qtot = None
        else:
            print(f"You didn't set the `qtot` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qtot={self.qtot}")
        
        if not hasattr(self, "qirreduced"):
            self.qirreduced = None
        else:
            print(f"You didn't set the `qirreduced` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qirreduced={self.qirreduced}")

        if not hasattr(self, "qinserted"):
            self.qinserted = None
        else:
            self.qinserted = int(self.qinserted)

        if not hasattr(self, "qirreduced_coords"):
            self.qirreduced_coords = None
        else:
            print(f"You didn't set the `qirreduced_coords` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0.")

        if not hasattr(self, "path_name_coords"):
            self.path_name_coords = None
        
        if self.mode == "matdyn":
            self.path_name_coords = self.get_hspp()

    def get_q_from_scfout(self, dir):
        if not os.path.exists(dir):
            raise FileExistsError ("scf.out didn't exist!")
        content = open(dir, "r").readlines()
        def find_k(item):
            if re.search(r"k\(\s*\d+\)\s*=\s*", item):
                return item

        result = filter(find_k, content)
        for res in result:
            ks  = re.findall(r"\-?\d+\.\d+", res.split(",")[0])
            self.q_coordinate_list.append(ks)
            wp = re.findall(r"\-?\d+\.\d+", res.split(",")[1])
            nqs = float(wp[0]) * self.qtot / 2
            self.q_weight_list.append(nqs)

        self.qtot           = self.qpoints[0] * self.qpoints[1] * self.qpoints[2] 
        self.q_non_irreducible_amount = len(self.q_coordinate_list)

        return self.qtot,    self.q_non_irreducible_amount, \
               self.q_coordinate_list, self.q_weight_list

    def get_qpoints(self, dyn0_path):
        '''
        input  : dir        dyn0文件的路径
                 
        return : self.qtot 总q点数
                 self.qirreduced 不可约q点数
                 self.qirreduced_coords 不可约q点坐标
        '''
        if not os.path.exists(dyn0_path):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dyn0_path, "r").readlines()
        # check qtot is right or not! 
        qpoints = list(map(int, content[0].strip("\n").split()))

        return qpoints

    def get_q_from_dyn0(self, dyn0_path):
        '''
        input  : dir        dyn0文件的路径
                 
        return : self.qtot 总q点数
                 self.qirreduced 不可约q点数
                 self.qirreduced_coords 不可约q点坐标
        '''
        if not os.path.exists(dyn0_path):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dyn0_path, "r").readlines()
        # check qtot is right or not! 
        qpoints = list(map(int, content[0].strip("\n").split()))
        qtot  = list(map(int, qpoints))
        if qtot != self.qpoints:
            print(f"qtot = {qtot}")
            print(f"self.qpoints = {self.qpoints}")
            raise ValueError ("q points set wrong")

        def find_q(item):
            if re.search(r"E[\+|\-]", item):
                return item

        qirreduced = int(content[1])
        _q_coordinate_list     = list(filter(find_q, content))
        qirreduced_coords = [q_string.strip("\n").split() for q_string in _q_coordinate_list]
        
        return qtot, qirreduced, qirreduced_coords

    def get_q_from_q2r(self, q2r_path):
        """
        input:  q2r_path   q2r.out文件的路径
        
        return:  self.qweights 不可约q点权重
        """
        def find_q_weight(item):
            if re.search(r"nqs\=", item):
                return item
        q2r_out        = open(q2r_path, "r").readlines()
        _q_weight_list = list(filter(find_q_weight, q2r_out))
        qweights  = [re.search(r"\d+", qwt).group() for qwt in _q_weight_list]
        
        return qweights

    def merge(self, dir):
        elph_dir_path = os.path.join(dir, "elph_dir")
        if not os.path.exists(elph_dir_path):
            os.makedirs(elph_dir_path)
        
        for i in range(self.qirreduced):
            src_elph   = os.path.join(dir, str(i+1), "elph_dir", "elph.inp_lambda.1")
            dst_elph   = os.path.join(elph_dir_path, "elph.inp_lambda."+str(i+1))
            shutil.copy(src_elph, dst_elph)
            print(f"elph.inp_.1 copy finished \n {dst_elph}")

            src_dyn    = os.path.join(dir, str(i+1), self.system_name+".dyn")
            dst_dyn    = os.path.join(dir,           self.system_name+".dyn"+str(i+1))
            shutil.copy(src_dyn, dst_dyn)
            print(f"{self.system_name}.dyn copy finished {dst_dyn}")

            for j in range(51, 61):
                src_a2Fq2r = os.path.join(dir, str(i+1), "elph_dir", "a2Fq2r."+str(j)+".1")
                dst_a2Fq2r = os.path.join(elph_dir_path,             "a2Fq2r."+str(j)+"."+str(i+1))
                shutil.copy(src_a2Fq2r, dst_a2Fq2r)
                print(f"a2Fq2r.{str(j)}.1 copy finished  {dst_dyn}")

    def get_hspp(self):
        """
        This method is to get high symmetry paths and points
        """ 
        lat             = self.ase_type.cell.get_bravais_lattice()
        path_name_string= lat.special_path

        _path_name_list = [[ p for p in pp] for pp in path_name_string.split(",")]
        print(f"the high symmetry points path is {_path_name_list}")
        high_symmetry_type = input("please input the mode you want: 'all_points', 'main_points', or Nothing to input\n\n")
        if not high_symmetry_type:
            high_symmetry_type = "all_points" # default
        if "," in path_name_string:
            if high_symmetry_type == "all_points":
                path_name_list = list(chain.from_iterable(_path_name_list))
                print(f"the choosed high symmetry points path is \n {path_name_list}")
            elif high_symmetry_type == "main_points":
                path_name_list = max(_path_name_list, key=len)
                print(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in path_name_string]
            print(f"the choosed high symmetry points path is \n {path_name_list}")

        special_points       = lat.get_special_points()
        path_coords          = [list(special_points[point_name]) for point_name in path_name_list]
        path_name_coords= list(zip(path_name_list, path_coords))

        print("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            print("{:<4} {:<10} {:<10} {:<10}".format(dirt[0], dirt[1], dirt[2], name))
        print("Print Cartesian coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            cart = np.dot(self.reciprocal_plattice, dirt)
            print("{:<10} {:<10} {:<10} {:<4} ".format(cart[0], cart[1], cart[2], name))

        return path_name_coords 

    def get_top_freq(self, dosfile):
        phonon_dos = open(dosfile, "r").readlines()
        frequents = []
        for id, line in enumerate(phonon_dos):
            if id != 0:
                linelist = line.split()
                frequents.append(float(linelist[0]))
        # convert cm-1 to Thz
        frequents = map(lambda x:x/33.3, frequents)
        top_freq  = ceil(max(frequents))

        return top_freq


class qedos_inputpara(qe_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict,
        ):
    
        super(qedos_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs
            )
        
        print('\n'); print("if you want to run `dos`, you had better set these values!")
        # 电子态密度设置
        if not hasattr(self, "DeltaE"):
            self.DeltaE = 0.01
            print("You didn't set `DeltaE` for eledos.in, the program will use default value: DeltaE=0.01")
        if not hasattr(self, "emin"):
            self.emin = -10
            print("You didn't set `emin`   for eledos.in, the program will use default value: emin=-10")
        if not hasattr(self, "emax"):
            self.emax = 30
            print("You didn't set `emax`   for eledos.in, the program will use default value: emax=30 ")

        # 声子态密度设置
        if hasattr(self, "qpoints"):
            _qpoints = self.qpoints.split()
            self.qpoints = list(map(int, _qpoints))
        else:
            dyn0_names = Path(self.work_underpressure).joinpath(f"{self.system_name}.dyn0")
            qpoints = self.get_qpoints(dyn0_path=dyn0_names)
            self.qpoints = [q*2 for q in qpoints]
            print("if you want to calculate phonodos, you had better set `qpoints`! ")
            print("----Its `qpoints` had better be set more densely than `qpoints` in `ph.in`")
            print("----For example, In ph.in, qpoints='8 8 8', then in phono_dos.in, qpoints='16 16 16' ")
            print(f"----You didn't set `qpoints`, the program will use `qpoints` in {self.system_name}.dyn0, and then multiply 2 for qpoints")

        if not hasattr(self, "ndos"):
            self.ndos = 500
            print("You didn't set `ndos`, the program will use default value: ndos=500")

    def get_qpoints(self, dyn0_path):
        '''
        input  : dir        dyn0文件的路径
                 
        return : self.qtot 总q点数
                 self.qirreduced 不可约q点数
                 self.qirreduced_coords 不可约q点坐标
        '''
        if not os.path.exists(dyn0_path):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dyn0_path, "r").readlines()
        # check qtot is right or not! 
        qpoints = list(map(int, content[0].strip("\n").split()))

        return qpoints


class qesc_inputpara(qephono_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict
    ):
        super(qesc_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs
            )
        print('\n'); print("if you want to run `superconduct`, you had better set these values!")
        # Mc-A-D and Eliashberg
        if not hasattr(self, "screen_constant"):
            print("You didn't set the `screen_constant` ! The program will use default value: screen_constant=0.13")
            self.screen_constant = 0.13

        # Mc-A-D
        q2r_names = list(Path(self.work_underpressure).glob("q2r.out"))
        if len(q2r_names)==1:
            q2r_path = q2r_names[0]
            self.qweights = self.get_q_from_q2r(q2r_path=q2r_path)
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            raise FileExistsError("No exist *.dyn0! ") 
        
        # Mc-A-D
        if not hasattr(self, "top_freq"):
            print(f"You didn't set the `top_freq`!         The program will read {self.system_name}_phono.dos file")
            phonodos_file = Path(self.work_underpressure).joinpath(self.system_name+"_phono.dos")
            if phonodos_file.exists():
                self.top_freq = self.get_top_freq(dosfile=phonodos_file)
        
        # Mc-A-D
        if not hasattr(self, "deguass"):
            self.deguass = 0.12
            print("You didn't set the `deguass` !         The program will use default value: deguass=0.12")

        # Mc-A-D
        if not hasattr(self, "smearing_method"):
            self.smearing_method = 1
            print("You didn't set the `smearing_method`!  The program will use default value: smearing_method=1")

        # Eliashberg
        print("If you use Eliashberg method, you have to specify the temperature_points !")
        print("If you use Eliashberg method, you may not specify the a2f_dos* !")
        print("If you use Eliashberg method, you may not specify the degauss_column* !")
        print("\tIf you set a2f_dos*,then you don't need set degauss_column !\n\tIf you set both, the program will run in the way of `degauss_column*`")
        if not hasattr(self, "temperature_points"):
            self.temperature_points = 5000
            print("You didn't set the `temperature_points`.The program will use default value: temperature_points=5000")

        # Eliashberg
        if not hasattr(self, "a2F_dos"):
            self.a2F_dos = None
            print("You didn't set the `a2F_dos`.          The program will use default value: a2F_dos=None")

        # Eliashberg
        if not hasattr(self, "degauss_column"):
            self.degauss_column = None
            print("You didn't set the `a2F_dos`.          The program will use default value: degauss_column=None")
        
        if hasattr(self, "Tc"):
            self.identify_converge_Tc()

    def identify_converge_Tc(self):
        '''
        根据四面体的非自洽方法计算得到的TDOS。
        然后取scf.out里面的费米能级。
        然后给TDOS中的横坐标能量减去一个费米能级得到最终可以出图的数据。(千万注意不能用.tdos中的费米能级, 它是qe做非自洽计算得到的。它并不准确。一定要取自洽计算得到的费米能级)

        无论是用Mc-A-D公式还是用Eliashberg公式计算超导, 都需要取一个合理的展宽对应的lambda对应的Tc
        那么如何得到这个合理的展宽呢???我们可以做非自洽计算得到TDOS, 找到费米能级处的总态密度, 然后对应找到该态密度在lamda.out文件中对应的费米能级处态密度。
        在lamda.out文件中对应的费米能级处态密度 还对应着一个lambda值, 这个lambda值在lamda.out中对应着一个Tc。这个Tc就是最终收敛的Tc.
        此时得到的Tc是通过Mc-A-D公式计算得到的。如果该Tc对应的lambda是一个大于1.5的值, 那么最好用eliashberg方程去求解Tc
            (注意: lamda.out中费米能级处的态密度是:  states/Ry/A3/Unit Cell/spin  )
            (注意: 但是qe用非自洽计算出的DOS的单位是: states/Ry/A3/Unit Cell       )

        那么如何得到 eliashberg方法中的Tc呢? 你需要手动选择一个合适的a2f.dos*. 
        如果你在计算声子的时候得到10个展宽, 那么在计算目录中计算phonodos时就会 得到10个a2f.dos*文件。
        通过前面计算非自洽的dos, 你可以在lamda.out中找到一个和非自洽dos的N(Ef)最接近的N(Ef), 这个N(Ef)对应的展宽就是一个合适的展宽。
        如果这个展宽是从大到小排列第7个展宽, 那么在选择a2f.dos*文 件时, 就选择a2f.dos7这个文件作为eliashberg的输入。

        但是现在有一个问题。为什么我在计算声子时，在ph.in文件中取了50个展宽，却只得到了10个a2f.dos文件

        '''
        # TODO 有待完善
        # tdos_file = Path(self.work_underpressure).joinpath(f"{self.system_name}.tdos")
        # if not tdos_file.exists():
            # print(f"{self.system_name}.tdos doesn't exist !!! The program will exit")
            # sys.exit(1)
        
        # 获得eliashberg方法计算得到的Tc
        print("运行可能报错, 建议仔细检查ELIASHBERG_GAP_T.OUT文件的第二列, 有些数据本来应该是0.1666489645E-265, 但是实际可能为0.1666489645-265, 导致numpy无法将其转化为数字。")
        eliashberg_gap_t_out = Path(self.work_underpressure).joinpath("ELIASHBERG_GAP_T.OUT")
        if not eliashberg_gap_t_out.exists():
            print(f"ELIASHBERG_GAP_T.OUT doesn't exist !!! The program will exit")
        gap_t   = np.loadtxt(eliashberg_gap_t_out)
        Tc      = gap_t[:, 0]
        gap     = gap_t[:, 1]
        dgap    = abs(np.gradient(gap, Tc))
        dgap_id = np.where(dgap < 1e-8)[0]
        gap_id  = np.where(gap  < 1e-8)[0]
        Tc_id   = [i for i in dgap_id if i in gap_id]
        print(f'Tc收敛温度为 {Tc[Tc_id[0]]} K (结果仅供参考, 请以Tc-GAP图像为准)')
        sys.exit(0)


class qeprepare_inputpara(qephono_inputpara):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict,
        ):

        super(qeprepare_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs
            )

        if hasattr(self, "mode"):
            self.mode = self.mode.split()