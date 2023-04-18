import os
import re
import sys
import shutil
import logging
from pathlib import Path
from itertools import chain
from math import ceil

import numpy as np
import pandas as pd

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
            print("Note: --------------------")
            print("    The program will set npool=1 in your submitjob scripts")
 
        if not hasattr(self, "queue"):
            self.queue = None
            print("Note: ----------------------")
            print("    You didn't specify queue, so the program will not submit the job in any way")

        print("Note: --------------------")
        print("if you want to run `relax` `scffit` `fit`, you had better set these values!")
        # &CONTROL
        if not hasattr(self, "forc_conv_thr"):
            self.forc_conv_thr = "1.0d-6"
            print("    You didn't set the `forc_conv_thr` ! The program will use default value: forc_conv_thr=1.0d-6")

        if not hasattr(self, "etot_conv_thr"):
            self.etot_conv_thr = "1.0d-7"
            print("    You didn't set the `etot_conv_thr` ! The program will use default value: etot_conv_thr=1.0d-7")

        # &SYSTEM
        if not hasattr(self, "occupations"):
            self.occupations = "smearing"
            print("    You didn't set the `occupations` !   The program will use default value: occupations=smearing")

        if not hasattr(self, "smearing"):
            self.smearing = "gauss"
            print("    You didn't set the `smearing` !      The program will use default value: smearing=gauss")

        # self.smearing = methfessel-paxton 做scffit 和 scf 时候用这个参数
        if not hasattr(self, "degauss"):
            self.degauss = "0.02"
            print("    You didn't set the `degauss` !       The program will use default value: degauss=0.02")

        if not hasattr(self, "ecutwfc"):
            self.ecutwfc = "60"
            print("    You didn't set the `ecutwfc` !       The program will use default value: ecutwfc=60")

        if not hasattr(self, "ecutrho"):
            self.ecutrho = "720"
            print("    You didn't set the `ecutrho` !       The program will use default value: ecutrho=720")

        if not hasattr(self, "lspinorb"):
            self.lspinorb = "false"
            print("    You didn't set the `lspinorb` !      The program will use default value: lspinorb=false")
        else:
            print("Please carefully check the bool value of `lspinorb` you just set. Its format must be `false` or `true` without capital")

        if self.lspinorb == "true":
            self.noncolin = "true"
            print("Because lspinorb = true, so the noncolin=true")
        elif not hasattr(self, "noncolin"):
            self.noncolin = "false"
            print("    You didn't set the `noncolin` !      The program will use default value: noncolin=false")
        else:
            print("Please carefully check the bool value of `noncolin` you just set. Its format must be `false` or `true` without capital")

        if not hasattr(self, "la2F"):
            self.la2F = "true"
            print("    You didn't set the `la2F` !          The program will use default value: la2F=true. ")
            print("                                         But in relax mode and scf mode, it doesn't exist ! It only exist in scffit mode")
        else:
            print("Please carefully check the bool value of `la2F` you just set. Its format must be `false` or `true` without capital")

        # &ELECTRONS
        if not hasattr(self, "diagonalization"):
            self.diagonalization = "david"
            print("    You didn't set the `diagonalization`! The program will use default value: diagonalization=david")
        
        if not hasattr(self, "conv_thr"):
            self.conv_thr = "1.0d-9"
            print("    You didn't set the `conv_thr` !      The program will use default value: conv_thr=1.0d-9")
            #  做结构弛豫 1.0-d8
            #  做scffit 和 scf 时候用1.0d-9

        if not hasattr(self, "mixing_beta"):
            self.mixing_beta = "0.7"
            print("    You didn't set the `mixing_beta` !   The program will use default value: mixing_beta=0.7")
            #  做scffit 和 scf 时候用0.8
        if not hasattr(self, "electron_maxstep"):
            self.electron_maxstep = "200"
            print("    You didn't set the `electron_maxstep`! The program will use default value: electron_maxstep=200")

        # &CELL
        if not hasattr(self, "press_conv_thr"):
            self.press_conv_thr = "0.01"
            print("    You didn't set the `press_conv_thr`! The program will use default value: press_conv_thr=0.01")

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
        
        print("Note: --------------------")
        print("if you want to run `phono`, you had better set these values!")
        dyn0_names = Path(self.work_path).joinpath(f"{self.system_name}.dyn0")
        if hasattr(self, "qpoints"):
            _qpoints = self.qpoints.split()
            self.qpoints = list(map(int, _qpoints))
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]
        elif dyn0_names.exists():
            print(f"    You didn't set the `qpoints` !        The program will qpoints in read {self.system_name}.dyn0 file")
            self.qpoints = self.get_qpoints(dyn0_path=dyn0_names)
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints] 
        else:
            self.qpoints = [4,4,4]
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]

        if not hasattr(self, "tr2_ph"):
            self.tr2_ph = "1.0d-16"
            print(f"    You didn't set the `tr2_ph` !          The program will use default value: tr2_ph=1.0d-16")

        if not hasattr(self, "electron_phonon"):
            self.electron_phonon="interpolated"
            print(f"    You didn't set the `electron_phonon` ! The program will use default value: electron_phonon=interpolated")

        if not hasattr(self, "el_ph_nsigma"):
            self.el_ph_nsigma = "10"
            print(f"You didn't set the `el_ph_nsigma` !    The program will use default value: el_ph_nsigma=10")
            
        if not hasattr(self, "el_ph_sigma"):
            self.el_ph_sigma = "0.005"
            print(f"    You didn't set the `el_ph_sigma` !     The program will use default value: el_ph_sigma=0.005")
        
        if not hasattr(self, "alpha_mix"):
            self.alpha_mix = "0.3"
            print(f"    You didn't set the `alpha_mix` !       The program will use default value: alpha_mix=0.3")
        
        if not hasattr(self, "dyn0_flag"):
            self.dyn0_flag = False
        
        else:
            self.dyn0_flag = eval(self.dyn0_flag)

        dyn0_names = Path(self.work_path).joinpath(f"{self.system_name}.dyn0")
        if dyn0_names.exists():
            self.qtot, self.qirreduced, self.qirreduced_coords= self.get_q_from_dyn0(dyn0_path=dyn0_names) 
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            print(f"{self.system_name}.dyn0 doesn't exist in {self.work_path}. The qtot, qirreduced, qirreduced_coords will not get values.")

        if not hasattr(self, "qtot"):
            self.qtot = None
        else:
            print(f"    You didn't set the `qtot` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qtot={self.qtot}")
        
        if not hasattr(self, "qirreduced"):
            self.qirreduced = None
        else:
            print(f"    You didn't set the `qirreduced` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qirreduced={self.qirreduced}")

        if not hasattr(self, "qinserted"):
            self.qinserted = None
        else:
            self.qinserted = int(self.qinserted)

        if not hasattr(self, "qirreduced_coords"):
            self.qirreduced_coords = None
        else:
            print(f"    You didn't set the `qirreduced_coords` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0.")

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
        _q_coordinate_list = list(filter(find_q, content))
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
        lat     = self.ase_type.cell.get_bravais_lattice()
        pstring = lat.special_path

        # 获得高对称点路径
        plist = [[ p for p in pp] for pp in pstring.split(",")]
        print(f"the high symmetry points path is {plist}")
        print(
            "please input the mode you want, just even input Number like 1 or 2\n",
            "'1  all_points'\n",
            "'2  main_points'\n",
            "Nothing to input, directly press ENTER, the default is main_points\n"
            )
        high_symmetry_type = input()
        if not high_symmetry_type:
            high_symmetry_type = "2" # default
        if "," in pstring:
            if high_symmetry_type == "1":
                path_name_list = list(chain.from_iterable(plist))
                print(f"the choosed high symmetry points path is \n {path_name_list}")
            elif high_symmetry_type == "2":
                path_name_list = max(plist, key=len)
                print(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in pstring]
            print(f"the choosed high symmetry points path is \n {path_name_list}")

        special_points   = lat.get_special_points()
        path_coords      = [list(special_points[point_name]) for point_name in path_name_list]
        path_name_coords = list(zip(path_name_list, path_coords))


        # 处理高对称点路径
        print("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            print("{:<10.6f} {:<10.6f} {:<10.6f} {:<4}".format(dirt[0], dirt[1], dirt[2], name))
        
        

        print("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
        for vector in self.reciprocal_plattice:
            print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))
                
        

        print("Print projected high symmetry path")
        projected_path_name_coords = [[path_name_coords[0][0], path_name_coords[0][1][0]]]
        total_dist = 0
        for idx in range(1, len(path_name_coords)):
            current_name   = path_name_coords[idx][0]
            current_coords = np.dot(self.reciprocal_plattice, path_name_coords[idx][1])
            last_coords    = np.dot(self.reciprocal_plattice, path_name_coords[idx-1][1])
            dist = np.linalg.norm(current_coords-last_coords, 2)
            total_dist += dist
            projected_path_name_coords.append([current_name, total_dist])
        string_names = '  '.join(coord[0] for coord in projected_path_name_coords)
        string_coord = '  '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
        print(string_names)
        print(string_coord)
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

    def check_convergence(self):
        """
        检查gaussian收敛性, 给出一个gaussian的索引, 获得收敛的gaussian
        """
        # 获得当前文件中关于gaussian值得设置
        # 检查 ph_no_split.in 文件是否存在，如果存在就，获取里面的展宽间距el_ph_sigma和展宽个数el_ph_nsigma
        ph_no_split_in_path = self.work_path.joinpath("ph_no_split.in")
        if ph_no_split_in_path.exists():
            greporder = f"grep el_ph_sigma {ph_no_split_in_path}"
            el_ph_sigma = float(re.search(r'\d+\.\d+', os.popen(greporder).read()).group())
            greporder = f"grep el_ph_nsigma {ph_no_split_in_path}"
            el_ph_nsigma = int(re.search(r'\d+', os.popen(greporder).read()).group())
            gaussian0 = [el_ph_sigma*(1+i) for i in range(0, el_ph_nsigma)]
        else:
            print("Note: --------------------")
            print("     The program didn't get `ph_no_split.in`. So you need input it by yourself.")
            el_ph_sigma = float(input("Input el_ph_sigma, it has to be a float number\n"))
            el_ph_nsigma= int(input("Input el_ph_nsigma, it has to be a integer number\n"))
            gaussian0 = [el_ph_sigma*(1+i) for i in range(0, el_ph_nsigma)]

        # 校验输入文件获得的gaussian和输出文件获得的gaussian是否一致
        greporder = f"grep Gaussian {self.work_path.joinpath('elph_dir', 'elph.inp_lambda.1')}" + "| awk '{print $3}'"
        gaussian1 = [float(num.strip("\n")) for num in os.popen(greporder).readlines()]
        flag = np.allclose(np.array(gaussian0), np.array(gaussian1), rtol=1e-6)
        if not flag:
            print("The Gaussion in ph_no_split.in file is different the Gaussion in `'elph_dir\elph.inp_lambda.1`")
            sys.exit(1)
        self.gaussian = gaussian1
        # 判断哪一个gaussian是收敛的，给出的是gaussian列表的索引值。
        shlorder = f"grep DOS {self.work_path.joinpath('elph_dir', 'elph.inp_lambda.1')}" + "| awk '{print $3}'"
        dos = os.popen(shlorder).readlines()
        dos = [float(i.split()[-1]) for i in dos]
        delta_dos = [abs(dos[i + 1] - dos[i]) for i in range(len(dos) - 1)]

        # 收敛Gaussian的索引
        idx = np.argmin(delta_dos)
        gauss = gaussian0[idx]
        print(f'Converged Index of gaussian={idx+1}, corresponding gaussian value={gaussian1[idx]}')

        return idx, gauss

    def get_gam_lines(self, gauss:float, q_number:int, freq_number:int):
        """
        从gam.lines文件中获得声子线宽
            gauss: 收敛的高斯展宽
            q_number: q点个数
            freq_number: 每个q点的振动模式数
        """

        # 从gam.lines文件中获得指定Gaussian展宽对应的所有的声子线宽
        gam_lines_path = self.work_path.joinpath("gam.lines")
        if not gam_lines_path.exists():
            print("Note: --------------------")
            print("    Sorry, The gam.lines doesn't exist !")
        tmp_gam1_path = self.work_path.joinpath("temp_gam1")
        tmp_gam2_path = self.work_path.joinpath("temp_gam2")
        print("Note: --------------------")
        print("    Two temporary files named `temp_gam1` ad `temp_gam2` will be created, which are based on gam.lines ")
        print("    They will be removed after obtaining the phonon-line-width of the corresponding Gaussian")
        os.system(f"""sed -n '/Broadening   {gauss}/,/Broadening/p' {gam_lines_path} > {tmp_gam1_path}""")
        # \w：匹配一个单词字符，包括字母、数字、下划线。
        # \+：匹配前面的字符或字符集至少一次或多次。
        # $：匹配输入的结尾。
        # \|：表示逻辑或，用于连接两个正则表达式。
        # Broadening：匹配字符串 "Broadening"。
        # .*：匹配任意数量的任意字符（除了换行符）。
        os.system(f"""sed 's/ \w\+$\|Broadening.*$//p' {tmp_gam1_path} > {tmp_gam2_path}""")
        phononwidth = []
        with open(tmp_gam2_path, "r") as temgam2:
            for line in temgam2.readlines():
                if line != "\n":
                    for pw in line.strip("\n").split():
                        phononwidth.append(pw)
        if len(phononwidth) != (q_number*freq_number):
            print("Note: --------------------")
            print("    The number of phonon-width is not equal to the number of `q_number*freq_number`")
            print(f"    phonon-width = {phononwidth}")
            print(f"    q_number * freq_number = {q_number} * {freq_number} = {q_number*freq_number}")
            sys.exit(1)
        # phononwidth 是一个一维数组，大小为 q点个数 * 每个q点的振动模式数
        phononwidth = np.array(phononwidth)
        phononwidth = phononwidth.reshape(q_number, freq_number)
        phononwidth = np.row_stack((phononwidth, [np.nan]*phononwidth.shape[1]))
        phononwidth = phononwidth.T.reshape(-1, 1) # 将其转化为一个 二维数组，只有一列。行数为：q点个数 * 每个q点的振动模式数
        os.system(f'rm -f {tmp_gam1_path} {tmp_gam2_path}')
        return phononwidth 
 
    def get_phono_freq(self):
        """获得可以在origin中作图的数据"""
        freq_gp_path = self.work_path.joinpath(self.system_name+".freq.gp")
        if not freq_gp_path.exists():
            print("Note: --------------------")
            print(f"    Sorry, {self.system_name}.freq.gp doesn't exist !")
        qpoints_freq = np.loadtxt(freq_gp_path)  # qpoints_freq 是一个 二维数组，第一列是q点坐标，其余列是振动频率
        q_number     = qpoints_freq.shape[0]     # q 点个数,          
        freq_number  = qpoints_freq.shape[1]-1   # 每个q点的振动模式数  -1是因为第一列是q点坐标，其余列是振动频率，所以要把第一列减去

        # 在最后增加一行nan用来 区别每一个q点对应的振动频率
        qpoints_freq = np.row_stack((qpoints_freq, [np.nan]*qpoints_freq.shape[1])) 
        
        qpoints = qpoints_freq[:, 0:1]  # 二维数组中取出一列qpoints_freq[:, 0]时，默认情况下会得到一个一维数组, 但是使用qpoints_freq[:, 0:1]就可以保持二维
        freqs   = qpoints_freq[:, 1:]   # 二维数组
        # 因为每个q点对应多个振动模式，需要将每个q点都复制f次，f是该q点对应的振动模式数
        # np.repeat 可以重复数组中的元素， q_points是一个q行1列的二维数组，q代表q点个数
        # 每个q点的重复次数为它的振动模式数，振动模式数即freq的列数freq.shape[1]，有freq.shape[1]列，就是有freq.shape[1]振动模式
        qpoints_new = np.repeat(qpoints, freqs.shape[1], axis=1) # 二维数组
        qpoints_new = qpoints_new.T.reshape(-1, 1)  # 将qpoints变成一个f*q行的二维数组
                                                    # 切记这里要对qpoints_new进行转置，因为reshape(-1,1)是取出来一整行，和下一整行进行拼接，不然就把一整行NaN拼接起来了。
        freqs_new   = freqs.T.reshape(-1, 1)        # 将freqs变成一个f*q行的二维数组
                                                    # 切记这里要对freqs进行转置，因为reshape(-1,1)是取出来一整行，和下一整行进行拼接，不然就把一整行NaN拼接起来了。
        qpoints_freq_new = np.hstack((qpoints_new, freqs_new)) # 沿着水平方向堆叠

        return qpoints_freq_new, q_number, freq_number

    def merge_qp_freq_width(self, qpoints_freqs, phononwidth):
        """
        合并q点位置，振动频率，声子线宽成一个3列n行的数组
        qpoints_freqs 是一个2列n行的数组
        phononwidth   是一个1列n行的数组
        """
        qpoints_freqs_phonowidth = np.hstack((qpoints_freqs, phononwidth))
        qp_freq_width = pd.DataFrame(qpoints_freqs_phonowidth)
        qp_freq_width.columns = ["qpoints", "freqs", "widths"]
        qp_freq_width.to_csv(
            self.work_path.joinpath("qp_freq_width.csv"),
            header=True,
            index=False,
            )
        

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
        
        print("Note: --------------------")
        print("if you want to run `dos`, you had better set these values!")
        # 电子态密度设置
        if not hasattr(self, "DeltaE"):
            self.DeltaE = 0.01
            print("    You didn't set `DeltaE` for eledos.in, the program will use default value: DeltaE=0.01")
        if not hasattr(self, "emin"):
            self.emin = -10
            print("    You didn't set `emin`   for eledos.in, the program will use default value: emin=-10")
        if not hasattr(self, "emax"):
            self.emax = 30
            print("    You didn't set `emax`   for eledos.in, the program will use default value: emax=30 ")

        # 声子态密度设置
        if hasattr(self, "qpoints"):
            _qpoints = self.qpoints.split()
            self.qpoints = list(map(int, _qpoints))
        else:
            dyn0_names = Path(self.work_path).joinpath(f"{self.system_name}.dyn0")
            qpoints = self.get_qpoints(dyn0_path=dyn0_names)
            self.qpoints = [q*2 for q in qpoints]
            print("Note: --------------------")
            print("if you want to calculate phonodos, you had better set `qpoints`! ")
            print("    Its `qpoints` had better be set more densely than `qpoints` in `ph.in`")
            print("    For example, In ph.in, qpoints='8 8 8', then in phono_dos.in, qpoints='16 16 16' ")
            print(f"    You didn't set `qpoints`, the program will use `qpoints` in {self.system_name}.dyn0, and then multiply 2 for qpoints")

        if not hasattr(self, "ndos"):
            self.ndos = 500
            print("    You didn't set `ndos`, the program will use default value: ndos=500")

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
    
    def get_phonodos(self):
        """
        返回值是一个dataframe类型，是一个二维列表，第一行是freq+元素名称, 第二行
        """
        # 获得phonodos计算的输出文件
        phonon_dos_path = self.work_path.joinpath(self.system_name+"_phono.dos")
        if not phonon_dos_path.exists():
            print("Note: --------------------")
            print(f"    Sorry, {self.system_name}_phono.dos doesn't exist !")
            sys.exit(1)

        # 获得元素的顺序 以及 每种元素的原子个数
        scffitin_path = self.work_path.joinpath("scffit.in")
        if not scffitin_path.exists():
            print("Note: --------------------")
            print(f"    Sorry, scffit.in doesn't exist !")
            sys.exit(1) 
        shlorder = "sed -n '/ATOMIC_POSITIONS (crystal)/,/K_POINTS {automatic}/ {//!p}' " + f"{scffitin_path}"
        elements_coords = os.popen(shlorder).readlines()
        elements = [ele.split()[0] for ele in elements_coords]

        phonondos = pd.read_table(
            phonon_dos_path,
            skiprows=1,  # skiprows=1：跳过文件的第一行，即不将其作为数据的一部分进行读取。
            header=None, # header=None：不将文件的第一行作为列名，而将其视为数据。
            sep='\s+'    # sep='\s+'：使用正则表达式 \s+ 作为列之间的分隔符，表示一个或多个空格字符。
            )
        columns = ['freq', 'tdos']+[ele for ele in elements]
        phonondos.columns = columns
        # 对相同的列名(相同的元素的pdos相加)进行相加,
        phonondos_sum = phonondos.groupby(phonondos.columns, axis=1).sum()
        # 将 'freq' 列移动到 DataFrame 的第一列
        freq_col = phonondos_sum.pop('freq')
        phonondos_sum.insert(0, 'freq', freq_col)
        # 将 'tdos' 列移动到 DataFrame 的第二列
        tdos_col = phonondos_sum.pop('tdos')
        phonondos_sum.insert(1, 'tdos', tdos_col)

        phonondos_sum.to_csv(
            self.work_path.joinpath("phonodos_proj2eles.csv"),
            header=True,   # 指定第一行为列名称
            index=False,   # 不输出列索引
            )

        return phonondos_sum


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
        print("Note: --------------------")
        print("if you want to run `superconduct`, you had better set these values!")
        # Mc-A-D and Eliashberg
        if not hasattr(self, "screen_constant"):
            print("    You didn't set the `screen_constant` ! The program will use default value: screen_constant=0.13")
            self.screen_constant = 0.13

        # Mc-A-D
        q2r_names = list(Path(self.work_path).glob("q2r.out"))
        if len(q2r_names)==1:
            q2r_path = q2r_names[0]
            self.qweights = self.get_q_from_q2r(q2r_path=q2r_path)
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            raise FileExistsError("No exist *.dyn0! ") 
        
        # Mc-A-D
        if not hasattr(self, "top_freq"):
            print(f"    You didn't set the `top_freq`!         The program will read {self.system_name}_phono.dos file")
            phonodos_file = Path(self.work_path).joinpath(self.system_name+"_phono.dos")
            if phonodos_file.exists():
                self.top_freq = self.get_top_freq(dosfile=phonodos_file)
        
        # Mc-A-D
        if not hasattr(self, "deguass"):
            self.deguass = 0.12
            print("    You didn't set the `deguass` !         The program will use default value: deguass=0.12")

        # Mc-A-D
        if not hasattr(self, "smearing_method"):
            self.smearing_method = 1
            print("    You didn't set the `smearing_method`!  The program will use default value: smearing_method=1")

        # Eliashberg
        print("If you use Eliashberg method, you have to specify the temperature_steps !")
        print("If you use Eliashberg method, you may not specify the a2f_dos* !")
        print("If you use Eliashberg method, you may not specify the degauss_column* !")
        print("\tIf you set a2f_dos*,then you don't need set degauss_column !\n\tIf you set both, the program will run in the way of `degauss_column*`")
        if not hasattr(self, "temperature_steps"):
            self.temperature_steps = 5000
            print("    You didn't set the `temperature_steps`.The program will use default value: temperature_steps=5000")

        # Eliashberg
        if not hasattr(self, "a2F_dos"):
            self.a2F_dos = None
            print("    You didn't set the `a2F_dos`.          The program will use default value: a2F_dos=None")

        # Eliashberg
        if not hasattr(self, "degauss_column"):
            self.degauss_column = None
            print("    You didn't set the `a2F_dos`.          The program will use default value: degauss_column=None")
        
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
        # tdos_file = Path(self.work_path).joinpath(f"{self.system_name}.tdos")
        # if not tdos_file.exists():
            # print(f"{self.system_name}.tdos doesn't exist !!! The program will exit")
            # sys.exit(1)
        
        # 获得eliashberg方法计算得到的Tc
        # 运行可能报错, 建议仔细检查ELIASHBERG_GAP_T.OUT文件的第二列, 有些数据本来应该是0.1666489645E-265, 但是实际可能为0.1666489645-265, 导致numpy无法将其转化为数字
        print("-------------------------------------WARING--------------------------------------")
        print("It is recommended to carefully check the second column of the ELIASHBERG_GAP_T.OUT file. ")
        print("Some data should be 0.1666489645E-265, but it may actually be 0.1666489645-265, causing numpy to fail to turn it into a number.")
        print("-------------------------------------WARING--------------------------------------")
        eliashberg_gap_t_out = Path(self.work_path).joinpath("ELIASHBERG_GAP_T.OUT")
        process_eliashberg_gap_t_out = Path(self.work_path).joinpath("process_ELIASHBERG_GAP_T.OUT")
        if not eliashberg_gap_t_out.exists():
            print(f"ELIASHBERG_GAP_T.OUT doesn't exist !!! The program will exit")
        os.system(f"sed '/E/!d' {eliashberg_gap_t_out} > {process_eliashberg_gap_t_out}") # 将不包含字符串E的行都删掉
        gap_t   = np.loadtxt(process_eliashberg_gap_t_out)
        Tc      = gap_t[:, 0]
        gap     = gap_t[:, 1]
        dgap    = abs(np.gradient(gap, Tc))
        dgap_id = np.where(dgap < 1e-8)[0]
        gap_id  = np.where(gap  < 1e-8)[0]
        Tc_id   = [i for i in dgap_id if i in gap_id]
        print(f'Tc = {Tc[Tc_id[0]]} K (The results are for reference only; please refer to the Tc-GAP image)')
        sys.exit(1)

    def calculate_lambda_personly(self, converged_gauss_index, gauss):
        """
        输入: 
            converged_degauss_index: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
            gauss: 对应索引的Gaussian值
        """
        from scipy.integrate import trapz

        alpha2F_dat_path = self.work_path.joinpath("alpha2F.dat")
        if not alpha2F_dat_path.exists():
            print("Note: --------------------")
            print("    `alpha2F.dat` file doesn't exist!!!")
            sys.exit(1)
        shlorder = f"sed -i 's/#//g' {alpha2F_dat_path}" # 这是因为第一行由于一个#, 删除它才不影响pandas读入 # E(THz)     0.005     0.010     0.015     0.020     0.025     0.030     0.035     0.040     0.045     0.050
        os.system(shlorder)
        alpha2F_dat = pd.read_table(
            alpha2F_dat_path,
            sep='\s+',
        )
        print("Note: --------------------")
        print(f"    Converged gauss index inputed = {converged_gauss_index+1}, its value = {gauss}")
        print(f"    So in alpha2F.dat, corresponding gauss value = {alpha2F_dat.columns[converged_gauss_index+1]} ")
        omegas  = np.array(alpha2F_dat.iloc[:, 0])
        alpha2f = np.array(alpha2F_dat.iloc[:, converged_gauss_index+1])

        # 定义被积函数的插值函数
        integrand_value = np.nan_to_num(alpha2f / omegas, nan=0) # # 计算 alpha2f / omegas，并将结果中的 NaN 替换为 0

        # 初始化结果列表
        integral_results = []
        # 遍历 omegas 数组，对每个积分上限进行数值积分
        for idx, omega in enumerate(omegas):
            # 使用 trapz 函数进行数值积分
            integral = trapz(y=integrand_value[:idx+1], x=omegas[:idx+1], dx=0.01)
            # 将积分结果添加到结果列表中
            integral_results.append(integral)
        print(integral_results)

    def obtain_lambda_omegalog_Tc_fromQE(self, converged_gauss_index, gauss):
        """
        输入: 
        converged_gauss_index: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
        """
        lambda_out_path = self.work_path.joinpath("lambda.out")
        if not lambda_out_path.exists():
            print("Note: --------------------")
            print("    `lambda.out` file doesn't exist!!!")
            sys.exit(1)
        shlorder = "sed -n '/^lambda/,$ {/^lambda/!p}'" + f"  {lambda_out_path}" # -n 选项表示只打印匹配的行
        content  = [line.strip('\n').split() for line in os.popen(shlorder).readlines()]
        degauss_lambda_omegalog_Tc = [[self.gaussian[idx], float(line[0]), float(line[1]), float(line[2])] for idx, line in enumerate(content)]
        print("Note: --------------------")
        print(f"    Converged Gaussian = {degauss_lambda_omegalog_Tc[converged_gauss_index][0]}")
        print(f"    Corresponding lambda = {degauss_lambda_omegalog_Tc[converged_gauss_index][1]}")
        print(f"    Corresponding omega_log = {degauss_lambda_omegalog_Tc[converged_gauss_index][2]}")
        print(f"    Corresponding Tc = {degauss_lambda_omegalog_Tc[converged_gauss_index][3]}") 


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