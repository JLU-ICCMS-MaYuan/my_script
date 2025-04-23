import os
import re
import sys
import logging
import shutil
from pathlib import Path
from itertools import chain
from math import ceil

import numpy as np
import pandas as pd

from qe.qe_base import qe_base

logger = logging.getLogger(__name__)

class qe_inputpara(qe_base):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict, # 这里非常重要, 因为qe_inputpara的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
        ):
        
        super(qe_inputpara, self).__init__(
            work_path,
            press,
            submit_job_system,
            input_file_path
        )
        
        logger.info('Run `inputpara`')

        for key, value in kwargs.items():
            setattr(self, key, value)
        
        if not hasattr(self, "mode"):
            logger.error("You must specify mode")
            sys.exit(1)

        if not hasattr(self, "execmd"):
            logger.warning("You must specify execute command, such as 'mpirun -np 48', 'bash', 'srun', 'srun --mpi=pmi2'")
            sys.exit(1)

        if not hasattr(self, "npool"):
            self.npool = 1
            logger.debug(f'npool = {self.npool}\n')
            logger.warning("You must specify npool, when you do SCTK and EPW calculations, you have to make sure npool = np")

        if not hasattr(self, "queue"):
            self.queue = None
            logger.debug(f'queue = {self.queue}\n')
        
        # &CONTROL
        if not hasattr(self, "forc_conv_thr"):
            self.forc_conv_thr = "1.0d-6"

        if not hasattr(self, "etot_conv_thr"):
            self.etot_conv_thr = "1.0d-7"

        # &SYSTEM
        if not hasattr(self, "occupations"):
            self.occupations = "smearing"

        if not hasattr(self, "smearing"):
            self.smearing = "gauss"

        if not hasattr(self, "nbnd"):
            self.nbnd = None

        # self.smearing = methfessel-paxton 做scffit 和 scf 时候用这个参数
        if not hasattr(self, "degauss"):
            self.degauss = "0.02"

        if not hasattr(self, "ecutwfc"):
            self.ecutwfc = "80"

        if not hasattr(self, "ecutrho"):
            self.ecutrho = "960"

        if not hasattr(self, "lspinorb"):
            self.lspinorb = "false"
        else:
            logger.debug(f'Please carefully check the bool value of `lspinorb` you just set. \nIts format must be `false` or `true` without capital\n')

        if self.lspinorb == "true":
            self.noncolin = "true"
            logger.debug(f'Because lspinorb = true, so the noncolin=true\n')
        elif not hasattr(self, "noncolin"):
            self.noncolin = "false"
            logger.debug(f"You didn't set the `noncolin` ! The program will use default value: noncolin=false\n")
        else:
            logger.debug(f"Please carefully check the bool value of `noncolin` you just set. Its format must be `false` or `true` without capital")

        if not hasattr(self, "la2F"):
            self.la2F = ".true."
            logger.debug("You didn't set the `la2F` ! The SCRIPTS4QE will use default value: la2F = .true. ")
            logger.debug("if you do nscf for SCTK, you have to set la2F = .true. But other nscf calculation, la2F = .false.")
            logger.debug("But in relax mode and scf mode, it doesn't exist ! It only exist in scffit mode")
        else:
            logger.debug("Please carefully check the bool value of `la2F` you just set. Its format must be `false` or `true` without capital")

        # &ELECTRONS
        if not hasattr(self, "diagonalization"):
            self.diagonalization = "david"
        
        if not hasattr(self, "conv_thr"):
            self.conv_thr = "1.0d-9"
            #  做结构弛豫 1.0-d8
            #  做scffit 和 scf 时候用1.0d-9

        if not hasattr(self, "mixing_beta"):
            self.mixing_beta = "0.7"
            #  做scffit 和 scf 时候用0.8
        
        if not hasattr(self, "electron_maxstep"):
            self.electron_maxstep = "200"
            logger.debug("You didn't set the `electron_maxstep`! The program will use default value: electron_maxstep=200")

        if not hasattr(self, "charge_density_dat"):
            self.charge_density_dat = "tmp/H3S1.save/charge-density.dat"

        if not hasattr(self, "data_file_schema_xml"):
            self.data_file_schema_xml = f"tmp/{self.system_name}.save/data-file-schema.xml"

        if not hasattr(self, "data_file_schema_xml"):
            self.data_file_schema_xml = f"tmp/{self.system_name}.save/data-file-schema.xml"

        if not hasattr(self, "paw_txt"):
            logger.warning("If you use PAW pseudopotential, you have to set `paw_txt`")
            self.paw_txt = f"tmp/{self.system_name}.save/paw.txt"

        # &CELL
        if not hasattr(self, "press_conv_thr"):
            self.press_conv_thr = "0.01"

        # &kpoints
        logger.debug("You have been confirmed that the kpoints_dense, kpoints_sparse, qpoints are all right and consistent with the symmetry")
        logger.debug("If you not make sure, you had better run 'vaspkit' or 'kmesh.py' to check it!!!")
        logger.debug('The order for vaspkit is:  (xxxx is corresponding to 1/LATTICE_PARA(i)*KPOINTS(i), (i) represents one axis, such as x_axis, y_axis, z_axis)')
        logger.debug(r'echo -e "1\n102\n2\nxxxx" | vaspkit')
        logger.debug('The order for kmesh.py is:  (xxxx is corresponding to the KSPACING in vasp)')
        logger.debug("kmesh.py xxxx")
        # time.sleep(3)
        if hasattr(self, "kpoints_dense"):
            _kpoints_dense = self.kpoints_dense.split()
            self.kpoints_dense = list(map(int, _kpoints_dense))
            logger.info('kpoints_dense={} by custom'.format(self.kpoints_dense))
        else:
            logger.warning('kpoints_dense={}'.format(None))


        if hasattr(self, "kpoints_sparse"):
            _kpoints_sparse = self.kpoints_sparse.split()
            self.kpoints_sparse = list(map(int, _kpoints_sparse))
            logger.info('kpoints_sparse={} by custom'.format(self.kpoints_sparse))
        else:
            logger.warning('kpoints_sparse={}'.format(None))

        if not hasattr(self, "wan"):
            self.wan = False
        else:
            self.wan = eval(self.wan)

        if not hasattr(self, "k_automatic"):
            self.k_automatic = True
            self.kpoints_coords = None
            self.totpts = 0
            logger.debug("Do not set k_automatic, default k_automatic = True")
        else:
            if self.k_automatic == "True":
                self.k_automatic = True; self.totpts = 0
                logger.debug("Set k_automatic = True by custom. Therefore, self.totpts = 0")
            else:
                self.k_automatic = False
                self.kpoints_coords, self.totpts = self.get_kmesh_justlike_kmesh_pl(kpoints=self.kpoints_dense)
                logger.debug("Set k_automatic = False by custom. Therefore, self.totpts was gotten by program")

        # 是否自动选择高对称路径
        if not hasattr(self, "autoselect"):
            self.autoselect = False
        else:
            self.autoselect = eval(self.autoselect)

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

    def get_hspp(self, autoselect=False):
        """
        This method is to get high symmetry paths and points
        """ 
        lat     = self.ase_type.cell.get_bravais_lattice()
        pstring = lat.special_path

        # 获得高对称点路径
        _plist  = [[ p for p in pp if not p.isdigit()] for pp in pstring.split(",")]
        logger.info(f"the high symmetry points path: \n{_plist}")

        print(
            "please input the mode you want, just even input Number like 1 or 2\n",
            "0:  all_points:\n",
            "1:  first_group_points\n",
            "2:  second_group_points\n",
            "n:  n^Th_group_points\n",
            "1 2: first_group_points and second_group_points"
            "...."
            "Nothing to input, directly press ENTER, the default is all_points\n"
            )
        if autoselect:
            high_symmetry_type = [0]
        else:
            high_symmetry_type = list(map(int, input().split())) #将输入的整数字符串按照空格划分成列表并分别转化为整数类型并再转化为列表
            if high_symmetry_type == []:
                high_symmetry_type = [0]
                logger.error("what you input is not an integer number, So use the `0:  all_points`")

        path_name_list = []
        if "," in pstring:
            if 0 in high_symmetry_type:
                path_name_list = list(chain.from_iterable(_plist))
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
            elif 0 not in high_symmetry_type:
                # path_name_list = max(_plist, key=len)
                for hst in high_symmetry_type:
                    path_name_list.extend(_plist[hst-1])
                logger.info(f"the choosed high symmetry points path is \n {path_name_list}")
        else:
            path_name_list = [ pp for pp in pstring]

        special_points   = lat.get_special_points()
        path_coords      = [list(special_points[point_name]) for point_name in path_name_list]
        path_name_coords = list(zip(path_name_list, path_coords))


        # 处理高对称点路径
        logger.info("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            print("{:<10.6f} {:<10.6f} {:<10.6f} {:<4}".format(dirt[0], dirt[1], dirt[2], name))
        
        

        logger.info("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
        for vector in self.reciprocal_plattice:
            print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))
                
        

        logger.info("Print projected high symmetry path")
        logger.info("倒格子的单位是 2pi/alat")
        #projected_path_name_coords = [[path_name_coords[0][0], path_name_coords[0][1][0]]]
        projected_path_name_coords = [[path_name_coords[0][0], 0]]
        total_dist = 0
        for idx in range(1, len(path_name_coords)):
            current_name   = path_name_coords[idx][0]
            # current_coords = np.dot(self.reciprocal_plattice, path_name_coords[idx][1])
            # last_coords    = np.dot(self.reciprocal_plattice, path_name_coords[idx-1][1])
            current_coords = np.dot(path_name_coords[idx][1],   self.reciprocal_plattice)
            last_coords    = np.dot(path_name_coords[idx-1][1], self.reciprocal_plattice)
            dist = np.linalg.norm(current_coords-last_coords, 2)
            total_dist += dist
            projected_path_name_coords.append([current_name, total_dist])
        string_names = ' '.join(coord[0] for coord in projected_path_name_coords)
        string_coord = ' '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
        logger.info(string_names)
        logger.info(string_coord)
        return path_name_coords 

    def get_kmesh_justlike_kmesh_pl(self, kpoints):
        """
        读取self.kpoints_dense参数, 将其传给n1, n2, n3, 再将n1, n2, n3转化为相应的倒空间的均匀网格点坐标
        """
        # 获取输入的 n1, n2, n3
        n1, n2, n3 = kpoints

        # 参数检查：确保 n1, n2, n3 都大于 0
        if n1 <= 0:
            print("n1 must be > 0")
            sys.exit()
        if n2 <= 0:
            print("n2 must be > 0")
            sys.exit()
        if n3 <= 0:
            print("n3 must be > 0")
            sys.exit()

        # 计算总的 k 点数量
        totpts = n1 * n2 * n3

        kpoints_coords = []
        if not self.wan: # 前三列写k点倒空间分数坐标，第四列写其权重
            print("K_POINTS crystal")
            print(totpts)
            for x in range(n1):
                for y in range(n2):
                    for z in range(n3):
                        # 格式化输出 k 点信息
                        # print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}{1/totpts:14.6e}")
                        kpoints_coords.append(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}{1/totpts:14.6e}")
        else:  # 只写前三列写k点倒空间分数坐标
            for x in range(n1):
                for y in range(n2):
                    for z in range(n3):
                        # 格式化输出 k 点信息（没有权重）
                        # print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}")
                        kpoints_coords.append(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}")
        return kpoints_coords, totpts

    
class qephono_inputpara(qe_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict, # 这里非常重要, 因为 qephono_inputpara 的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
        ):
        super(qephono_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs, # 这里非常重要, 因为qephono_inputpara的__init__需要读入kwargs
            )
        
        logger.info("Run `phono`")
        
        dyn0_names = Path(self.work_path).joinpath(f"{self.system_name}.dyn0")
        if hasattr(self, "qpoints"):
            _qpoints = self.qpoints.split()
            self.qpoints = list(map(int, _qpoints))
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]
            logger.info('qpoints={}, kpoints_sparse={}, kpoints_dense={} by custom'.format(self.qpoints, self.kpoints_sparse, self.kpoints_dense))
        elif dyn0_names.exists():
            self.qpoints = self.get_qpoints(dyn0_path=dyn0_names)
            self.kpoints_sparse= [kp*2 for kp in self.qpoints]
            self.kpoints_dense = [kp*4 for kp in self.qpoints]
            logger.info('qpoints={}, kpoints_sparse={}, kpoints_dense={} by dyn0'.format(self.qpoints, self.kpoints_sparse, self.kpoints_dense))
        else:
            logger.warning('qpoints={}, kpoints_sparse={}, kpoints_dense={}'.format(None, None, None))

        if not hasattr(self, "EPC_flag"):
            self.EPC_flag = True
        else:
            self.EPC_flag = eval(self.EPC_flag)

        if not hasattr(self, "search_sym"):
            self.search_sym = ".true." # 当计算SCTK时一般会设置为false
            
        # 针对SCTK计算的参数----------------------------
        self.SCTK_flag = False
        logger.debug("If you wanna use SCTK phonon and EPC, you have to use class qesctk_inputpara")
        
        if not hasattr(self, "tr2_ph"):
            self.tr2_ph = "1.0d-14"

        if not hasattr(self, "electron_phonon"):
            self.electron_phonon="interpolated"

        if not hasattr(self, "el_ph_nsigma"):
            self.el_ph_nsigma = "10"
            
        if not hasattr(self, "el_ph_sigma"):
            self.el_ph_sigma = "0.005"
        
        if not hasattr(self, "alpha_mix"):
            self.alpha_mix = "0.3"
        
        if not hasattr(self, "trans"):
            self.trans = ".true."
        
        if not hasattr(self, "dyn0_flag"):
            self.dyn0_flag = False
        else:
            self.dyn0_flag = eval(self.dyn0_flag)
        
        dyn0_names = Path(self.work_path).joinpath(f"{self.system_name}.dyn0")
        if dyn0_names.exists():
            self.qtot, self.qirreduced, self.qirreduced_coords= self.get_q_from_dyn0(dyn0_path=dyn0_names) 
            logger.debug("qtot from No.1 line in {}.dyn0 = {}".format(self.system_name, self.qtot))
            logger.debug("qirreduced from No.2 line in  {}.dyn0 = {}".format(self.system_name, self.qirreduced))
            logger.debug("qirreduced_coords from the rest lines in {}.dyn0, the number = {}".format(self.system_name, len(self.qirreduced_coords)))
            # time.sleep(3)
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            logger.error(f"{self.system_name}.dyn0 doesn't exist in {self.work_path}. The qtot, qirreduced, qirreduced_coords will not get values.")

        if not hasattr(self, "qtot"):
            self.qtot = None
        else:
            logger.debug(f"You didn't set the `qtot` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qtot={self.qtot}")
        
        if not hasattr(self, "qirreduced"):
            self.qirreduced = None
        else:
            logger.debug(f"You didn't set the `qirreduced` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0. qirreduced={self.qirreduced}")

        if not hasattr(self, "qinserted"):
            self.qinserted = 50
        else:
            self.qinserted = int(self.qinserted)

        if not hasattr(self, "qirreduced_coords"):
            self.qirreduced_coords = None
        else:
            logger.debug(f"You didn't set the `qirreduced_coords` ! Of course you don't need to set it! The program will get from {self.system_name}.dyn0.")
        
        # phonodos计算
        if not hasattr(self, "ndos"):
            self.ndos = 500

        # phonoband计算
        if self.mode == "matdyn" or self.mode == "processphono":
            self.path_name_coords = self.get_hspp(autoselect=self.autoselect)
        else:
            self.path_name_coords = None

        # 这一部分是关于如何获得收敛的gaussid和gauss
        if self.mode == "Tc" or self.mode == "phonobandwidthsdata":
            if hasattr(self, "gaussid") and hasattr(self, "gauss"):
                self.gaussid = int(self.gaussid)
                self.gauss = float(self.gauss)
                logger.info(f'Converged gaussid = {self.gaussid+1}, corresponding gaussian value={self.gauss}')
            elif hasattr(self, "efermi_dos"):
                self.efermi_dos = float(self.efermi_dos)
                self.gaussid, self.gauss = self.check_convergence(efermi_dos=self.efermi_dos)
                logger.info(f'Converged gaussid = {self.gaussid+1}, corresponding gaussian value={self.gauss}')
            else:
                self.efermi_dos = None
                self.gaussid, self.gauss = self.check_convergence(efermi_dos=self.efermi_dos)
                logger.info(f'Converged gaussid = {self.gaussid+1}, corresponding gaussian value={self.gauss}')

    def get_qpoints(self, dyn0_path):
        '''
        input  : directory        dyn0文件的路径
                 
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
        input  : directory        dyn0文件的路径
                 
        return : self.qtot 总q点数
                 self.qirreduced 不可约q点数
                 self.qirreduced_coords 不可约q点坐标
        '''
        if not os.path.exists(dyn0_path):
            raise FileExistsError ("dyn0 file doesn't exist!")
        content = open(dyn0_path, "r").readlines()
        # check qtot is right or not! 
        qpoints   = list(map(int, content[0].strip("\n").split()))
        qtot_list = list(map(int, qpoints))
        qtot      = qtot_list[0] * qtot_list[1] * qtot_list[2]

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

    def merge(self, directory):
        elph_dir_path = os.path.join(directory, "elph_dir")
        if not os.path.exists(elph_dir_path):
            os.makedirs(elph_dir_path)
        
        for i in range(self.qirreduced):

            src_split_in   = os.path.join(directory, str(i+1), "split_ph.in")
            dst_split_in   = os.path.join(directory, "split_ph.in."+str(i+1))
            src_split_out  = os.path.join(directory, str(i+1), "split_ph.out")
            dst_split_out  = os.path.join(directory, "split_ph.out."+str(i+1))
            if os.path.exists(src_split_in) and os.path.exists(src_split_out):
                shutil.copy(src_split_in, dst_split_in)
                logger.debug(f"split_ph.in copy finished \n {src_split_in}")
                shutil.copy(src_split_out, dst_split_out)
                logger.debug(f"split_ph.out copy finished \n {src_split_in}")
            else:
                logger.error(f"In {str(i+1)}, split_ph.in or split_ph.out doesn't exist ! Exit the program!")
                sys.exit(1)

            src_elph_NoNum = os.path.join(directory, str(i+1), "elph_dir", "elph.inp_lambda.1")
            src_elph_Num   = os.path.join(directory, str(i+1), "elph_dir", "elph.inp_lambda."+str(i+1))
            dst_elph = os.path.join(elph_dir_path, "elph.inp_lambda."+str(i+1))
            if os.path.exists(src_elph_NoNum):
                src_elph = src_elph_NoNum
                shutil.copy(src_elph, dst_elph)
                logger.debug(f"elph.inp_.1 copy finished \n {dst_elph}")
            elif  os.path.exists(src_elph_Num):
                src_elph = src_elph_Num
                shutil.copy(src_elph, dst_elph)
                logger.debug(f"elph.inp_.1 copy finished \n {dst_elph}")
            else:
                logger.error(f"In {str(i+1)}, elph.inp_lambda-file doesn't exist ! Exit the program!")


            


            src_dyn_NoNum= os.path.join(directory, str(i+1), self.system_name+".dyn")
            src_dyn_Num  = os.path.join(directory, str(i+1), self.system_name+f".dyn{str(i+1)}")
            if os.path.exists(src_dyn_NoNum):
                src_dyn = src_dyn_NoNum
            elif  os.path.exists(src_dyn_Num):
                src_dyn = src_dyn_Num
            else:
                logger.error(f"In {str(i+1)}, dyn-file doesn't exist ! Exit the program!")
                sys.exit(1)
            dst_dyn    = os.path.join(directory, self.system_name+".dyn"+str(i+1))
            shutil.copy(src_dyn, dst_dyn)
            logger.debug(f"{self.system_name}.dyn copy finished {dst_dyn}")

            for j in range(51, 51+int(self.el_ph_nsigma)):
                src_a2Fq2r_NoNum = os.path.join(directory, str(i+1), "elph_dir", "a2Fq2r."+str(j)+".1")
                src_a2Fq2r_Num   = os.path.join(directory, str(i+1), "elph_dir", "a2Fq2r."+str(j)+"."+str(i+1))
                dst_a2Fq2r = os.path.join(elph_dir_path, "a2Fq2r."+str(j)+"."+str(i+1))
                if os.path.exists(src_a2Fq2r_NoNum):
                    src_a2Fq2r = src_a2Fq2r_NoNum
                    shutil.copy(src_a2Fq2r, dst_a2Fq2r)
                    logger.debug(f"a2Fq2r.{str(j)}.1 copy finished  {dst_a2Fq2r}")
                elif  os.path.exists(src_a2Fq2r_Num):
                    src_a2Fq2r = src_a2Fq2r_Num
                    shutil.copy(src_a2Fq2r, dst_a2Fq2r)
                    logger.debug(f"a2Fq2r.{str(j)}.1 copy finished  {dst_a2Fq2r}")
                else:
                    logger.error(f"{src_a2Fq2r_NoNum} and {src_a2Fq2r_Num} doesn't exist ! Exit the program!")


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

    def check_convergence(self, efermi_dos=None):
        """
        检查gaussian收敛性, 给出一个gaussian的索引, 获得收敛的gaussian
        """
        # 获得当前文件中关于gaussian值得设置
        # 检查 ph_no_split.in 文件是否存在，如果存在就，获取里面的展宽间距el_ph_sigma和展宽个数el_ph_nsigma
        ph_no_split_in_path = self.work_path.joinpath("ph_no_split.in")
        if ph_no_split_in_path.exists():
            # greporder = f"grep el_ph_sigma {ph_no_split_in_path}"
            # el_ph_sigma = float(re.search(r'\d+\.\d+', os.popen(greporder).read()).group())
            # greporder = f"grep el_ph_nsigma {ph_no_split_in_path}"
            # el_ph_nsigma = int(re.search(r'\d+', os.popen(greporder).read()).group())
            el_ph_sigma = float(self.el_ph_sigma)
            el_ph_nsigma = int(self.el_ph_nsigma)
            gaussian0 = [el_ph_sigma*(1+i) for i in range(0, el_ph_nsigma)]
        else:
            logger.info("The program didn't get `ph_no_split.in`. So you need input it by yourself.")
            el_ph_sigma = float(input("Input el_ph_sigma, it has to be a float number\n"))
            el_ph_nsigma= int(input("Input el_ph_nsigma, it has to be a integer number\n"))
            gaussian0 = [el_ph_sigma*(1+i) for i in range(0, el_ph_nsigma)]

        # 校验输入文件获得的gaussian和输出文件获得的gaussian是否一致
        if not self.work_path.joinpath('elph_dir', 'elph.inp_lambda.1').exists():
            logger.error("    elph_dir/elph.inp_lambda.1 doesn't exist !!! The program will exit!!!")
            sys.exit(1)
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
        if not dos:
            logger.error("grep DOS elph_dir/elph.inp_lambda.1 failes !!! The progam will exit !!! Maybe calculation in elph get something wrong!!!")
            sys.exit(1)
        
        dos = [float(i.strip('\n')) for i in dos]
        if efermi_dos is None:
            logger.debug("You didn't set efermi_dos, so the program will use its own method to get converged gauss value.")
            delta_dos = [np.abs(dos[i + 1] - dos[i]) for i in range(len(dos) - 1)]
            idx = np.argmin(delta_dos) # 收敛Gaussian的索引
            gauss = gaussian0[idx]
        else:
            delta_dos = [np.abs(idos-float(efermi_dos)) for idos in dos]
            idx = np.argmin(delta_dos) # 收敛Gaussian的索引
            gauss = gaussian0[idx]
        # time.sleep(3)
        return idx+1, gauss  # 因为python都是从0开始

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
            logger.error("Sorry, The gam.lines doesn't exist !")
        tmp_gam1_path = self.work_path.joinpath("temp_gam1")
        tmp_gam2_path = self.work_path.joinpath("temp_gam2")
        logger.debug("Two temporary files named `temp_gam1` ad `temp_gam2` will be created, which are based on gam.lines ")
        logger.debug("They will be removed after obtaining the phonon-line-width of the corresponding Gaussian")
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
            logger.info("The number of phonon-width is not equal to the number of `q_number*freq_number`")
            logger.debug(f"phonon-width = {phononwidth}")
            logger.info(f"q_number * freq_number = {q_number} * {freq_number} = {q_number*freq_number}")
            sys.exit(1)
        # phononwidth 是一个一维数组，大小为 q点个数 * 每个q点的振动模式数
        phononwidth = np.array(phononwidth)
        phononwidth = phononwidth.reshape(q_number, freq_number)
        phononwidth = np.row_stack((phononwidth, [np.nan]*phononwidth.shape[1]))
        phononwidth = phononwidth.T.reshape(-1, 1) # 将其转化为一个 二维数组，只有一列。行数为：q点个数 * 每个q点的振动模式数
        os.system(f'rm -f {tmp_gam1_path} {tmp_gam2_path}')
        return phononwidth 
 
    def read_hspp_in_matdyn(self):

        if not self.work_path.joinpath("matdyn.in").exists():
            logger.warning(f"There is no matdyn.in in {self.work_path}")
            return None
        
        # 都当前目录下的matdyn.in
        matdyn_in_file = self.work_path.joinpath("matdyn.in")
        with open(matdyn_in_file, 'r') as f:
            lines = f.readlines()
            # 获取\所在的行号
            for i, line in enumerate(lines, start=1):
                if '/' in line:
                    slash_line_number = i
                    break
            else:
                logger.warning("There is non '/' in matdyn.in")
                return None

        # 获得高对称点路径
        path_name_coords = []
        for line in lines[slash_line_number+1:]:
            coords = list(map(float, line.split()[0:3]))
            name = line.split('!')[-1].strip('\n')
            path_name_coords.append([name, coords])

        # 获得倒格子
        logger.debug("Print projected high symmetry path")
        logger.debug("倒格子的单位是 2pi/alat")
        logger.debug("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
        for vector in self.reciprocal_plattice:
            print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))
        
        # 处理高对称点路径
        logger.debug("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            print("{:<10.6f} {:<10.6f} {:<10.6f} {:<4}".format(dirt[0], dirt[1], dirt[2], name))
        
        
        logger.info("Print projected high symmetry path")
        logger.info("倒格子的单位是 2pi/alat")
        #projected_path_name_coords = [[path_name_coords[0][0], path_name_coords[0][1][0]]]
        projected_path_name_coords = [[path_name_coords[0][0], 0]]
        total_dist = 0
        for idx in range(1, len(path_name_coords)):
            current_name   = path_name_coords[idx][0]
            # current_coords = np.dot(self.reciprocal_plattice, path_name_coords[idx][1])
            # last_coords    = np.dot(self.reciprocal_plattice, path_name_coords[idx-1][1])
            current_coords = np.dot(path_name_coords[idx][1],   self.reciprocal_plattice)
            last_coords    = np.dot(path_name_coords[idx-1][1], self.reciprocal_plattice)
            dist = np.linalg.norm(current_coords-last_coords, 2)
            total_dist += dist
            projected_path_name_coords.append([current_name, total_dist])
        string_names = ' '.join(coord[0] for coord in projected_path_name_coords)
        string_coord = ' '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
        logger.info(string_names)
        logger.info(string_coord)

    def get_phono_freq(self):
        """获得可以在origin中作图的数据"""
        freq_gp_path = self.work_path.joinpath(self.system_name+".freq.gp")
        if not freq_gp_path.exists():
            logger.warning(f"Sorry, {self.system_name}.freq.gp doesn't exist !")
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
        qp_freq_width.columns = ["Q-Path", "Frequency(cm\+(-1))", "Widths"]
        qp_freq_width.to_csv(
            self.work_path.joinpath("qp_freq_width.csv"),
            header=True,
            index=False,
            )
        
    def get_phonodos(self):
        """
        返回值是一个dataframe类型，是一个二维列表，第一行是freq+元素名称, 第二行
        """
        # 获得phonodos计算的输出文件
        phonon_dos_path = self.work_path.joinpath(self.system_name+"_phono.dos")
        if not phonon_dos_path.exists():
            logger.error(f"Sorry, {self.system_name}_phono.dos doesn't exist !")
            sys.exit(1)

        # 获得元素的顺序 以及 每种元素的原子个数
        logger.debug("Get the order of elements by the order of atomic coords in `scffit.in`")
        logger.debug("The order of elements in `***_phono.dos` is consistent with the order of `ATOMIC_SPECIES` in `scffit.in`")
        logger.debug("So, for guaranteeing the consistency of order,  the order of atomic coords in `scffit.in` has to be same as the order of `ATOMIC_SPECIES` in `scffit.in`")
        scffitin_path = self.work_path.joinpath("scffit.in")
        if not scffitin_path.exists():
            logger.error(f"Sorry, scffit.in doesn't exist !")
            sys.exit(1) 
        shlorder = "sed -n '/ATOMIC_POSITIONS/,/K_POINTS/ {//!p}' " + f"{scffitin_path}"        
        elements_coords = os.popen(shlorder).readlines()
        elements = [ele.split()[0] for ele in elements_coords if ele.strip()] # ele.strip() 作用是去掉 ATOMIC_POSITIONS 和 K_POINTS 之间的空行
        phonondos = pd.read_table(
            phonon_dos_path,
            skiprows=1,  # skiprows=1：跳过文件的第一行，即不将其作为数据的一部分进行读取。
            header=None, # header=None：不将文件的第一行作为列名，而将其视为数据。
            sep='\s+'    # sep='\s+'：使用正则表达式 \s+ 作为列之间的分隔符，表示一个或多个空格字符。
            )
        columns = ['freq', 'tdos']+[ele for ele in elements]
        phonondos.columns = columns
        # 对相同的列名(相同的元素的pdos相加)进行相加,
        # phonondos_sum = phonondos.groupby(phonondos.columns, axis=1).sum()
        phonondos_sum = phonondos.T.groupby(phonondos.columns).sum().T
        # 将 'freq' 列移动到 DataFrame 的第一列
        freq_col = phonondos_sum.pop('freq')
        phonondos_sum.insert(0, 'freq', freq_col)
        # 将 'tdos' 列移动到 DataFrame 的第二列
        tdos_col = phonondos_sum.pop('tdos')
        phonondos_sum.insert(1, 'tdos', tdos_col)

        phonondos_sum.to_csv(
            self.work_path.joinpath("phdos_proj2eles.csv"),
            header=True,   # 指定第一行为列名称
            index=False,   # 不输出列索引
            )

        return phonondos_sum

    def get_gibbs_from_phtdos(self):

        freq_phtdos_phpdos = self.get_phonodos()
        freq_phtdos_phpdos = freq_phtdos_phpdos[freq_phtdos_phpdos["freq"]>0]  # 只保留正频率
       
        freqs = freq_phtdos_phpdos['freq'].values # 讲pandas中提取出的freqs转化为numpy类型  # freqs is circular frequency
        freqs = freqs * 2.997924E+10 # convert unit cm-1 to Hz

        # 这种归一化方式非常有趣，可以避免乘以频率步长d_freq
        logger.debug("检验是否将声子态密度归一化到3N") 
        sumdos = freq_phtdos_phpdos['tdos'].sum()
        freq_phtdos_phpdos['tdos'] = freq_phtdos_phpdos['tdos']/sumdos*3*self.all_atoms_quantity
        # 定义积分步长
        d_freq = np.round(0.333564E-10*(freqs[-1] - freqs[0])/(len(freqs) -1), decimals=6)  # conserve unit is cm-1
        logger.debug("归一化至3N前:", d_freq*sumdos)
        logger.debug("归一化至3N后:", freq_phtdos_phpdos['tdos'].sum())

        tdos   = freq_phtdos_phpdos['tdos'].values  # conserve unit is states/cm-1 = states·cm
        temperature = np.array([i for i in range(0,5100,100)])

        # 定义计算单个声子的吉布斯自由能gibbs_i 
        def gibbs_i(freqs, T, dos):
            """\int{ G_i * g(w) } = \int{ (zpe_i + temperature_effect_i ) * g(w) dw}"""
            # h_bar = 6.582120E-16 # unit is eV·s
            h     = 4.13566770E-15 #  unit is eV·s 普朗克常数
            k_B   = 8.61734315E-5  # unit is eV/K
            if T == 0:
                gibbs_i = (0.5 * h * freqs) * dos
            else:
                gibbs_i = (0.5 * h * freqs + k_B * T * np.log(1-np.exp(-(h*freqs)/(k_B*T)))) * dos
                # gibbs_i = k_B*T * np.log(2*np.sinh(h*freqs/(2*k_B*T))) * dos
            return gibbs_i

        T_gibbs = []
        for T in temperature:
            gibbs = gibbs_i(freqs, T, tdos) #* d_freq # 这里因为前面特殊的归一化方式避免了乘以步长的问题
            gibbs = np.sum(gibbs)
            T_gibbs.append([T, gibbs])

        with open(self.work_path.joinpath("thermodynamics_from_phtdos.csv"), "w") as f:
            f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
            for T, gibbs in T_gibbs:
                f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/self.all_atoms_quantity))

    def get_gibbs_from_freq(self):
        """
        从声子总态密度获得声子频率然后获得自由能
            声子态密度的单位是: states/cm-1, 通过对声子态密度积分就可以得到总的振动模式数: 3N, 其中N是胞内总原子数
        """
        temperature = np.array([i for i in range(0,5100,100)])

        def gibbs_q(freqs, T, weight):
            """\int{ G_i * g(w) } = \int{ (zpe_i + temperature_effect_i ) * g(w) dw}"""
            h_bar = 6.582119514E-16 # unit is eV·s
            k_B   = 8.617343E-5  # unit is eV/K
            freqs = freqs[ freqs > 0 ] * 2*np.pi # 圆频率
            # print(len(freqs)); print(weight)
            # input(freqs)
            if T == 0:
                _gibbs_q = (1/2 * h_bar * freqs)
            else:
                _gibbs_q = (1/2 * h_bar * freqs + k_B * T * np.log(1-np.exp(-(h_bar*freqs)/(k_B*T))))
            # print(_gibbs_q); input()
            _gibbs_q = np.sum(_gibbs_q) * weight
            # print(_gibbs_q); input()
            return _gibbs_q
        dyn_paths = list(
            filter(lambda x: 'dyn' in x \
                        and 'dyn0' not in x \
                        and 'dyna2F' not in x 
                        and 'matdyn.modes' not in x \
                        and 'thermodynamics_from_dyn1.csv' not in x \
                        and "thermodynamics_from_dyns.csv" not in x \
                        and "thermodynamics_from_phtdos.csv" not in x, 
                        os.listdir(self.work_path))
                )        
        full_dyns_paths = sorted([os.path.join(self.work_path, dyn_path) for dyn_path in dyn_paths])
        if not full_dyns_paths:
            logger.warning(f'Unable to detect *dyn* files in {self.__path}. The program will exit')
            sys.exit(1)
        
        dyns  = [Dyn(path) for path in full_dyns_paths]

        # 从所有的q点获得频率计算吉布斯自由能
        T_gibbs_dyns = []
        for T in temperature:
            gibbs = 0
            for dyn in dyns:
                gibbs_Q = gibbs_q(dyn.freqs*1.0E+12, T, dyn.weight) # dyn.freqs 的单位是THz
                gibbs += gibbs_Q
            gibbs = gibbs/self.qtot # 除以总的q点数
            T_gibbs_dyns.append([T, gibbs])
        with open(self.work_path.joinpath("thermodynamics_from_dyns.csv"), "w") as f:
            f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
            for T, gibbs in T_gibbs_dyns:
                f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/self.all_atoms_quantity))

        # 从gamma点获得频率计算吉布斯自由能
        T_gibbs_dyn1 = []
        for T in temperature:
            gibbs = 0
            gibbs_Q = gibbs_q(dyns[0].freqs*1.0E+12, T, dyns[0].weight) # dyn.freqs 的单位是THz
            gibbs += gibbs_Q
            T_gibbs_dyn1.append([T, gibbs])
        with open(self.work_path.joinpath("thermodynamics_from_dyn1.csv"), "w") as f:
            f.write("{:>5},{:>12},{:>12}\n".format("T", "gibbs(eV)", "gibbs(eV/atom)"))
            for T, gibbs in T_gibbs_dyn1:
                f.write("{:>5},{:>12.8f},{:>12.8f}\n".format(T, gibbs, gibbs/self.all_atoms_quantity))


class Dyn(object):
    """
    Parses a single *dyn*.elph* file,
    returning all it contains:
    q-point, lambdas, gammas and squared frequencies
    """
    lines = list()
    q_point = tuple()
    weight = float()


    def __init__(self, path):
        with open(path) as read_obj:
            lines = read_obj.readlines()
            read_obj.close()
        self.lines = lines
        # print(f'q = ({", ".join(["%.3f" % round(_q, 3) for _q in self.q_point])}) '
            #   f'with number of q in the star {int(self.weight)}')
        
    @property
    def q_point(self):
        q_idx = int()
        for idx, line in enumerate(self.lines):
            if 'q = (' in line:
                q_idx = idx
                break
        self._q_point = tuple(float(x) for x in self.lines[q_idx].split()[-4:-1])
        return self._q_point
    
    @property
    def weight(self):
        self._weight = 0
        for line in self.lines:
            if 'Diagonalizing the dynamical matrix' in line: # 这一行前面都是等价的q点坐标，这一行后面是针对这些等价的坐标选取一个进行动力学矩阵对角化
                break
            if 'q = (' in line:
                self._weight = self._weight + 1
        # self.weight = float(self.lines[2].split()[1])
        return self._weight

    @property
    def freqs(self):
        self._freqs = []
        for line in self.lines:
            if 'freq (' in line:
                freq = float(line.split()[4])
                self._freqs.append(freq)
        self._freqs = np.array(self._freqs)
        return self._freqs


class qeeletron_inputpara(qe_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict,
        ):
        super(qeeletron_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs  # 这里非常重要, 因为 qeeletron_inputpara 的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
            )

        logger.info("Run `eletron`")
        # 电子态密度设置
        if not hasattr(self, "DeltaE"):
            self.DeltaE = 0.01
            logger.debug("You didn't set `DeltaE` for eledos.in, the program will use default value: DeltaE=0.01")
        if not hasattr(self, "emin"):
            self.emin = -10
            logger.debug("You didn't set `emin`   for eledos.in, the program will use default value: emin=-10")
        if not hasattr(self, "emax"):
            self.emax = 30
            logger.debug("You didn't set `emax`   for eledos.in, the program will use default value: emax=30. ")
            logger.debug("The reason is to avoid E-fermi bigger than then max-energy. ")
            logger.debug("For example, E-fermi=16eV, the upper energy you set is 10eV, then you can't get positive energys.")
        if not hasattr(self, "ngauss"):
            self.ngauss = 0
            logger.debug("You didn't set `emax`   for eledos.in, the program will use default value: ngauss=0 ")
        if not hasattr(self, "pdosdegauss"):
            self.pdosdegauss = 0
            logger.debug("You didn't set `emax`   for eledos.in, the program will use default value: pdosdegauss=0 ")
        # 电子能带计算
        if self.mode == "eleband" or self.mode == "eleproperties":
            self.path_name_coords = self.get_hspp(autoselect=self.autoselect)
        else:
            self.path_name_coords = None

        if not hasattr(self, "kinserted"):
            self.kinserted = 100
        else:
            self.kinserted = int(self.kinserted)
            
        if not hasattr(self, "nbnd"):
            logger.debug("You didn't set nbnd, the default value nbnd=100")
            self.nbnd = 100

        if not hasattr(self, "lsym"):
            logger.debug("You didn't set lsym, the default value lsym=.false.")
            self.lsym = 'false'


class qesc_inputpara(qephono_inputpara):

    def __init__(
        self, 
        work_path: str, 
        press: int, 
        submit_job_system: str, 
        input_file_path: str, 
        **kwargs: dict, 
    ):
        super(qesc_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs, # 这里非常重要, 因为 qesc_inputpara 的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
            )
        logger.info("run `superconduct`")
        # Mc-A-D and Eliashberg
        self.screen_constant = [0.10, 0.13, 0.16]

        # Mc-A-D
        q2r_names = list(Path(self.work_path).glob("q2r.out"))
        if len(q2r_names)==1:
            q2r_path = q2r_names[0]
            self.qweights = self.get_q_from_q2r(q2r_path=q2r_path)
            logger.debug("The number of qweights = {}".format(len(self.qweights)))
            # 获得 self.qtot, self.qirreduced, self.qirreduced_coords, self.qweights 
        else:
            logger.error("No exist q2r.out ! The program will exit !!!")
            sys.exit(1)
        
        # Mc
        phonodos_file = Path(self.work_path).joinpath(self.system_name+"_phono.dos")
        if phonodos_file.exists():
            self.top_freq = self.get_top_freq(dosfile=phonodos_file)
        else:
            logger.error(f"There is no {self.system_name}+'_phono.dos'")
            sys.exit(1)
        
        # Mc
        if not hasattr(self, "broaden"):
            self.broaden = 0.5
            logger.debug("You didn't set the `broaden` ! ")
            logger.debug("The program will use default value: broaden=0.5")
            logger.debug("***** Don't change easily the parameter of broaden, 0.5 is good enough !!!! ***** ")

        # Mc
        if not hasattr(self, "smearing_method"):
            self.smearing_method = 1
            logger.debug("You didn't set the `smearing_method`!  The program will use default value: smearing_method=1")

        # Eliashberg
        logger.debug("If you use Eliashberg method, you have to specify the temperature_steps, a2fdos*(alpha2fdat*)!")
        logger.debug("If you set both a2fdos* and alpha2fdat*, the program will run in the way of `a2fdos*`")
        if not hasattr(self, "temperature_steps"):
            self.temperature_steps = 500
            logger.debug("You didn't set the `temperature_steps`.The program will use default value: temperature_steps=5000")

        if hasattr(self, "a2fdos") and eval(self.a2fdos) == True:
            self.a2fdos = True
            self.alpha2fdat = False
            logger.debug("The omegas-alpha2F values will be getted from a2F.dos{}".format("xxx"))
            logger.debug("The shell order used is:")
            logger.debug("sed '1,5d' a2F.dos%s | sed '/lambda/d' | awk '{print $1/2, $2}' > ALPHA2F.OUT"%("xxx"))
        elif hasattr(self, "alpha2fdat") and eval(self.alpha2fdat) == True:
            self.a2fdos = False
            self.alpha2fdat = True
            logger.debug("The omegas-alpha2F values will be getted from alpha2F.dat{}".format("xxx"))
            logger.debug("The shell order used is:")
            logger.debug("sed '1,1d' alpha2F.dat | awk '{print $1/6579.684, $%s}' > ALPHA2F.OUT "%("xxx"))
        else:
            self.a2fdos = True
            self.alpha2fdat = False
            logger.debug("You didn't specify the either `a2fdos` or `alpha2fdat`")
            logger.debug("The omegas-alpha2F values willbe getted from a2F.dos")
            logger.debug("Because alpha2f.dat usually get something wrong")
            logger.debug("The shell order used is:")
            logger.debug("sed '1,5d' a2F.dos%s | sed '/lambda/d' | awk '{print $1/2, $2}' > ALPHA2F.OUT"%("xxx"))


    def get_a2Fdos_data(self, gauss_idx):
        a2Fdos_path = self.work_path.joinpath("a2F.dos"+str(gauss_idx))
        if not a2Fdos_path.exists():
            logger.error(f"    `a2F.dos+{str(gauss_idx)}` file doesn't exist!!!")
            sys.exit(1)
        # 这里简单说明一下a2F.dos*的文件结构
        # 频率             总a2F   剩下的列是N*3列，是N个原子，每个原子3个振动模式。
        # 特别注意第一列频率的单位是：里德伯Ry
        a2Fdos_data = os.popen(f"sed '1,5d' {a2Fdos_path} | sed '/lambda/d' ").readlines()
        a2Fdos_data = np.array([[float(x) for x in line.strip('\n').split()] for line in a2Fdos_data])

        return a2Fdos_data

    def get_alpha2Fdat_data(self, gauss_idx, gauss):
        alpha2F_dat_path = self.work_path.joinpath("alpha2F.dat")
        if not alpha2F_dat_path.exists():
            logger.error("`alpha2F.dat` file doesn't exist!!!")
            sys.exit(1)

        # 判断el_ph_nsigma是否大于10, 如果大于10,  alpha2F.dat文件折叠出2行来记录数据，巨恶心
        interval = ceil(int(self.el_ph_nsigma)/10)
        #  alpha2F.dat文件折叠出2行
        if interval == 2: 
            # 读取alpha2F.dat内容
            with open(alpha2F_dat_path, 'r') as file:
                lines = file.readlines()
            # 每n行合并成一行, n是el_ph_nsigma/10向上取整的结果
            merged_lines = []
            for i in range(0, len(lines), interval):
                merged_line = lines[i].strip() + " " + lines[i + 1].strip()
                merged_lines.append(merged_line)
            # 将合并的内容写入新的文件或直接读取
            merged_alpha2F_dat_path = self.work_path.joinpath('merged_alpha2F.dat')
            with open(merged_alpha2F_dat_path, 'w') as file:
                file.write("\n".join(merged_lines))
            shlorder = f"sed -i 's/#//g' {merged_alpha2F_dat_path}"; os.system(shlorder)
            shlorder = f"sed -i 's/#//g' {alpha2F_dat_path}"; os.system(shlorder)
            alpha2F_data = pd.read_table(
            merged_alpha2F_dat_path,
            sep='\s+',
            )
        elif interval == 1:
            shlorder = f"sed -i 's/#//g' {alpha2F_dat_path}" # 这是因为第一行由于一个#, 删除它才不影响pandas读入 # E(THz)     0.005     0.010     0.015     0.020     0.025     0.030     0.035     0.040     0.045     0.050
            os.system(shlorder)
            alpha2F_data = pd.read_table(
            alpha2F_dat_path,
            sep='\s+',
            )

        logger.debug(f"Converged gauss index inputed = {gauss_idx}, its value = {gauss}")
        logger.debug(f"So in alpha2F.dat, corresponding gauss value = {alpha2F_data.columns[gauss_idx]} ") # alpha2F_data.columns[gauss_idx]获取列的名称, 第0列是E(Thz)列, 第1列是0.005列. 由于我们输入的gauss_idx是从1开始的, 所以要获得正确的列不需要对gauss_idx做任何处理
        
        return alpha2F_data

    def get_lambda_from_alpha2f_single_broadening(self, alpha2Fdat_data, gauss_idx, gauss):
        """
        输入: 
            gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
            gauss: 对应索引的Gaussian值
        """
        
        # 计算lambda
        frequency  = np.array(alpha2Fdat_data.iloc[:, 0]); frequency[frequency == 0.0] = 0.00001
        a2F = np.array(alpha2Fdat_data.iloc[:, gauss_idx]) # alpha2Fdat_data 第0列是E(Thz)列, 第1列是0.005列. 由于我们输入的gauss_idx是从1开始的, 所以要获得正确的列相应的a2F不需要对gauss_idx做任何处理
        lambda_value = np.trapz(2 * a2F / frequency, frequency)

        # 将 frequency 和 a2F 作为列来创建 DataFrame
        w_alpha2f = pd.DataFrame({
            'omegas(Thz)': frequency,
            'alpha2f': a2F
        })
        w_alpha2f.to_csv(
            self.work_path.joinpath("w_alpha2f_from_alpha2Fdat.csv"),
            header=True,
            index=False,
        )
    
        return lambda_value

    def get_lambda_from_a2fdos_single_broadening(self, a2Fdos_data, gauss_idx):
        """
        输入: 
            gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
            gauss: 对应索引的Gaussian值
        """

        # 计算lambda
        frequency = a2Fdos_data[:,0]
        a2F = a2Fdos_data[:,1]
        lambda_value = np.trapz(2 * a2F / frequency, frequency)

        # 将 frequency 和 a2F 作为列来创建 DataFrame
        w_alpha2f = pd.DataFrame({
            'omegas(Rydberg)': frequency,
            'alpha2f': a2F
        })
        w_alpha2f.to_csv(
            self.work_path.joinpath("w_alpha2f_from_a2Fdos"+str(gauss_idx)+".csv"),
            header=True,
            index=False,
        )
        return lambda_value
    
    def get_wlog_from_a2fdos_single_broadening(self, a2Fdos_data, Lambda_bya2Fdos):
        """
        输入: 
            gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
            gauss: 对应索引的Gaussian值
        """

        # 计算wlog
        Ry2K = 157887.51240116
        frequency = a2Fdos_data[:,0]
        a2F = a2Fdos_data[:,1]
        w_log = Ry2K * np.exp(2/Lambda_bya2Fdos * np.trapz(a2F / frequency * np.log(frequency), frequency))
        return w_log

    def get_w2_from_a2fdos_single_broadening(self, a2Fdos_data, Lambda_bya2Fdos):
        """
        输入: 
            gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
            gauss: 对应索引的Gaussian值
        """
        # 计算w2
        frequency = a2Fdos_data[:,0]
        a2F = a2Fdos_data[:,1]
        w2  = np.sqrt(2/Lambda_bya2Fdos * np.trapz(frequency*a2F, frequency))
        return w2

    def getTc_by_eliashberg(self):
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
        logger.debug("It is recommended to carefully check the second column of the ELIASHBERG_GAP_T.OUT file. ")
        logger.debug("Some data should be 0.1666489645E-265, but it may actually be 0.1666489645-265, causing numpy to fail to turn it into a number.")
        eliashberg_gap_t_out = Path(self.work_path).joinpath("ELIASHBERG_GAP_T.OUT")
        process_eliashberg_gap_t_out = Path(self.work_path).joinpath("process_ELIASHBERG_GAP_T.OUT")
        if not eliashberg_gap_t_out.exists():
            logger.error(f"ELIASHBERG_GAP_T.OUT doesn't exist !!! The program will exit")
            sys.exit(1)
        os.system(f"sed '/E/!d' {eliashberg_gap_t_out} > {process_eliashberg_gap_t_out}") # 将不包含字符串E的行都删掉
        gap_t   = np.loadtxt(process_eliashberg_gap_t_out)
        Tc      = gap_t[:, 0]
        gap     = gap_t[:, 1]
        dgap    = abs(np.gradient(gap, Tc))
        dgap_id = np.where(dgap < 1e-8)[0]
        gap_id  = np.where(gap  < 1e-8)[0]
        Tc_id   = [i for i in dgap_id if i in gap_id]
        try:
            Tc = Tc[Tc_id[0]]
            return Tc
        except:
            logger.error("    Maybe the imaginary frequency of phono leads to NAN in ELIASHBERG_GAP_T.OUT. So The program will exit.")
            sys.exit(1)

    def getTc_McM_byqe(self, gauss_idx):
        """
        输入: 
        gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
        """
        lambda_out_path = self.work_path.joinpath("lambda.out")
        if not lambda_out_path.exists():
            logger.error("`lambda.out` file doesn't exist!!!")
            sys.exit(1)
        shlorder = "sed -n '/^lambda/,$ {/^lambda/!p}'" + f"  {lambda_out_path}" # -n 选项表示只打印匹配的行
        content  = [line.strip('\n').split() for line in os.popen(shlorder).readlines()]
        lambda_omegalog_Tc = [[float(line[0]), float(line[1]), float(line[2])] for idx, line in enumerate(content)]

        Lambda = lambda_omegalog_Tc[gauss_idx-1][0]
        omega_log = lambda_omegalog_Tc[gauss_idx-1][1]
        tc = lambda_omegalog_Tc[gauss_idx-1][2]
        return Lambda, omega_log, tc

    def getTc_by_McAD_from_a2fdos_single_broadening(
            self, 
            Lambda,
            wlog,
            w2,
            screen_constant):
        """
        输入: 
        gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
        """


        f1 = np.cbrt(1 + np.power(Lambda / (2.46 * (1 + 3.8 * screen_constant)), 1.5) )
        # f2 = 1 + (Lambda**2 * (w2 / wlog - 1)) / (Lambda**2 + 3.312 * ((1 + 6.3 * screen_constant) * w2 / wlog)**2) # it has issue in Allen-Dynes article.
        f2 = 1 - (Lambda**2 * (1-w2/wlog)) / (Lambda**2 + 3.312*(1+6.3*screen_constant)**2)
        Tc_McM = wlog/1.2 * np.exp( (-1.04*(1+Lambda)) / (Lambda-screen_constant*(1+0.62*Lambda)) )
        Tc_AD  = f1*f2*Tc_McM
        return Tc_McM, Tc_AD


class qebatch_inputpara(qephono_inputpara):

    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict, 
        ):

        super(qebatch_inputpara, self).__init__(
            work_path, 
            press, 
            submit_job_system, 
            input_file_path, 
            **kwargs, # 这里非常重要, 因为 qeprepare_inputpara 的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
            )
        logger.info("run `batch`")


class qesctk_inputpara(qephono_inputpara):
    
    def __init__(
        self,
        work_path: str,
        press: int,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict,
    ):
        super(qesctk_inputpara, self).__init__(
            work_path,
            press,
            submit_job_system,
            input_file_path,
            **kwargs # 这里非常重要, 因为qephono_inputpara的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
        )
        logger.info("run `SCTK`")
            

        self.SCTK_flag = True

        self.kpoints_coords_for_Twin, self.totpts_for_Twin = self.get_kmesh_justlike_twingrid_x(kpoints=self.qpoints)
        
        if not hasattr(self, "elph_nbnd_min"):
            self.elph_nbnd_min = 0
        
        if not hasattr(self, "elph_nbnd_max"):
            self.elph_nbnd_max = 0

        if not hasattr(self, "nbnd"):
            self.nbnd = None
            logger.warning("nbnd in `nscf.in` and `twin.in` have to be larger than that in `scf.in`")

        if not hasattr(self, "nci"):
            self.nci = 5
        
        if not hasattr(self, "laddxc"):
            self.laddxc = 0
            
        if not hasattr(self, "lsf"):
            self.lsf = 1
            
        if not hasattr(self, "ecutfock"):
            self.ecutfock = 80
            
        if not hasattr(self, "temp"):
            self.temp = -0.1
        
        if not hasattr(self, "fbee"):
            self.fbee = 1
            
        if not hasattr(self, "lbee"):
            self.lbee = 15
            logger.warning("Usually, in `kel.in`, lbee=15, in `lambda_mu_k, scdft_tc, deltaf, qpdos`, lbee=10.")
        
        if not hasattr(self, "xic"):
            self.xic = -1.0
        
        if not hasattr(self, "nmf"):
            self.nmf = 10
        
        if not hasattr(self, "nx"):
            self.nx = 100
        
        if not hasattr(self, "ne"):
            self.ne = 100

        if not hasattr(self, "emin"):
            self.emin = 1.0e-7
        
        if not hasattr(self, "emax"):
            self.emax = 0.7
            logger.warning("Usually, in `kel.in`, emax=5.0, in `lambda_mu_k, scdft_tc, deltaf, qpdos`, emax=0.7")

        if not hasattr(self, "electron_maxstep"):
            self.electron_maxstep = 100
        
        if not hasattr(self, "conv_thr"):
            self.conv_thr = 1.0e-15
        
        if not hasattr(self, "spin_fluc"):
            self.spin_fluc = ".true."

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
    
    def get_kmesh_justlike_twingrid_x(self, kpoints):
        # Print the header
        print("K_POINTS crystal")
        
        n1, n2, n3 = kpoints
        totpts = n1 * n2 * n3 * 2

        kpoints_coords = []
        # First loop (original points)
        for i1 in range(n1):
            kv1 = i1 / n1
            if (i1 * 2) >= n1:
                kv1 -= 1.0
            for i2 in range(n2):
                kv2 = i2 / n2
                if (i2 * 2) >= n2:
                    kv2 -= 1.0
                for i3 in range(n3):
                    kv3 = i3 / n3
                    if (i3 * 2) >= n3:
                        kv3 -= 1.0
                    # print(f"{kv1:20.15f} {kv2:20.15f} {kv3:20.15f}   1.0")
                    kpoints_coords.append(f"{kv1:20.15f} {kv2:20.15f} {kv3:20.15f}   1.0")

        # Second loop (special k-points)
        for i1 in range(n1):
            kv1 = (2 * i1 + 1) / (2 * n1)
            if (i1 * 2 + 1) >= n1:
                kv1 -= 1.0
            for i2 in range(n2):
                kv2 = (2 * i2 + 1) / (2 * n2)
                if (i2 * 2 + 1) >= n2:
                    kv2 -= 1.0
                for i3 in range(n3):
                    kv3 = (2 * i3 + 1) / (2 * n3)
                    if (i3 * 2 + 1) >= n3:
                        kv3 -= 1.0
                    # print(f"{kv1:20.15f} {kv2:20.15f} {kv3:20.15f}   1.0")
                    kpoints_coords.append(f"{kv1:20.15f} {kv2:20.15f} {kv3:20.15f}   1.0")
                    
        return kpoints_coords, totpts