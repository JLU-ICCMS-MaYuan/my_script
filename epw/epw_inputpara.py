from epw.epw_base import epw_base

import numpy as np

import logging
import sys
import os
from pathlib import Path

logger = logging.getLogger(__name__)

class epw_inputpara(epw_base):

    def __init__(
        self,
        work_path: str,
        submit_job_system: str,
        input_file_path: str,
        **kwargs: dict, # 这里非常重要, 因为 epw_inputpara 的__init__需要读入kwargs, 其它需要继承这个类的子类也需要保有这个参数kwargs
        ):
        super(epw_inputpara, self).__init__(
            work_path,
            submit_job_system,
            input_file_path,
            # 这里没有写 **kwargs, 是因为我继承了_base这个类, 这个类没有kwargs参数, 
            # 因为 epw_inputpara 的__init__需要初始化很多参数与的不同, 所以直接继承_base, 而没有继承_inputpara
        )
        logger.info("run `EPW`")
        
        for key, value in kwargs.items():
            setattr(self, key, value)

        if not hasattr(self, "execmd"):
            logger.debug("You must specify execute command, such as 'mpirun -np 48', 'bash', 'srun', 'srun --mpi=pmi2'")
            sys.exit(1)

        if not hasattr(self, "npool"):
            self.npool = 1
            logger.debug(f'npool = {self.npool}\n')
        logger.debug("You must guarantee that `npool` equal `cores, namly, -np`")

        if not hasattr(self, "queue"):
            self.queue = None
            logger.debug(f'queue = {self.queue}\n')
        
        if not hasattr(self, "mode"):
            logger.error("You must specify mode, such as 'epw_eband', 'epw_phono', 'epw_elph', 'epw_aniso_sc'")
            sys.exit(1)
            
        self.path_name_coords_for_EPW = self.get_hspp_for_EPW()

        if not hasattr(self, "dvscf_dir"):
            logger.error("You must specify dvscf_dir")
            sys.exit(1)
        else:
            if not Path.cwd().joinpath(self.dvscf_dir):
                logger.error(f"Specified dvscf_dir: {self.dvscf_dir} doesn't exist ")
                sys.exit(1)
            else:
                logger.info(f"Specified dvscf_dir: {self.dvscf_dir} exist.  ")
        logger.info(f"You have to write it as `save/`, not `save`. If you input reletive path, then you have to input `../save/`, not `../save`.")
        
        if not hasattr(self, "etf_mem"):
            self.etf_mem = 1
        
        if not hasattr(self, "lifc"):
            self.lifc = ".true."
        
        if not hasattr(self, "wannierize"):
            self.wannierize = '.false.'
            logger.debug(f'wannierize = {self.wannierize}\n')
        
        if not hasattr(self, "nbndsub"):
            logger.error(f"You must specify nbndsub!")
            sys.exit(1)
        else:
            logger.debug(f'nbndsub = {self.nbndsub}\n')

        if not hasattr(self, "bands_skipped"):
            logger.warning(f"You had better specify bands_skipped! For example, bands_skipped='1:10'")
            self.bands_skipped = None
        else:
            logger.debug(f'bands_skipped = {self.bands_skipped}\n')
        
        if not hasattr(self, "num_iter"):
            self.num_iter = 500

        if not hasattr(self, "dis_num_iter"):
            self.dis_num_iter = 200

        if not hasattr(self, "dis_froz_min"):
            logger.error(f"You must specify dis_froz_min !")
            sys.exit(1)
        else:
            logger.debug(f'dis_froz_min = {self.dis_froz_min}\n')

        if not hasattr(self, "dis_froz_max"):
            logger.error(f"You must specify dis_froz_max !")
            sys.exit(1)
        else:
            logger.debug(f'dis_froz_max = {self.dis_froz_max}\n')

        if not hasattr(self, "dis_win_max"):
            logger.error(f"You must specify dis_win_max !")
            sys.exit(1)
        else:
            logger.debug(f'dis_win_max = {self.dis_win_max}\n')

        if not hasattr(self, "proj"):
            logger.error(f"You must specify proj! Just like: proj='Nb:d  H:s'  (Attention: No space before or after the equal sign)")
            sys.exit(1)
        else:
            self.proj = self.proj.split()
            for idx, pj in enumerate(self.proj):
                logger.debug(f'proj({idx+1}) = {pj}\n')
        
        if not hasattr(self, "filkf") or not hasattr(self, "filqf"):
            logger.info("You didn't specify filkf or filqf, so the program will get it from modified_{}_band.kpt that is obtained from {}_band.kpt".format(self.system_name, self.system_name)) 
            self.filkf = "modified_{}_band.kpt".format(self.system_name)
            self.filqf = "modified_{}_band.kpt".format(self.system_name)
            old_band_kpt = self.work_path.joinpath("{}_band.kpt".format(self.system_name))
            new_band_kpt = self.work_path.joinpath("modified_{}_band.kpt".format(self.system_name))
            os.system("sed '1s/$/    crystal/' {} >  {}".format(old_band_kpt, new_band_kpt))

        if not hasattr(self, "asr_typ"):
            self.asr_typ = "simple"
            logger.debug(f'asr_typ = {self.asr_typ}\n')
        
        if not hasattr(self, "fsthick"):
            self.fsthick = 10
            logger.debug(f'fsthick = {self.fsthick}\n')
        
        if not hasattr(self, "degaussw"):
            self.degaussw = 0.1
            logger.debug(f'degaussw = {self.degaussw}\n')
        
        if not hasattr(self, "degaussq"):
            self.degaussq = 0.05
            logger.debug(f'degaussq = {self.degaussq}\n')
        
        if not hasattr(self, "nk"):
            logger.warning(f"You must specify nk! Just like: nk='4 4 4'  (Attention: No space before or after the equal sign)")
            self.nk = [4,4,4]
        else:
            self.nk = self.check_kgrid_qgird(self.nk)
            logger.debug(f'nk = {self.nk}\n')
        
        if not hasattr(self, "nq"):
            logger.warning(f"You must specify nk! Just like: nk='4 4 4'  (Attention: No space before or after the equal sign)")
            self.nq = [4,4,4]
        else:
            self.nq = self.check_kgrid_qgird(self.nq)
            logger.debug(f'nq = {self.nq}\n')
            
        if not hasattr(self, "nkf"):
            logger.warning(f"You must specify nfk! Just like: nfk='4 4 4'  (Attention: No space before or after the equal sign)")
            self.nkf = [1,1,1]
        else:
            self.nkf = self.check_kgrid_qgird(self.nkf)
            logger.debug(f'nkf = {self.nkf}\n')
            
        if not hasattr(self, "nqf"):
            logger.warning(f"You must specify nqf! Just like: nqf='4 4 4'  (Attention: No space before or after the equal sign)")
            self.nqf = [1,1,1]
        else:
            self.nqf = self.check_kgrid_qgird(self.nqf)
            logger.debug(f'nqf = {self.nqf}\n')
    
        if not hasattr(self, "wscut"):
            self.wscut = 3.0
            logger.debug(f'wscut = {self.wscut}\n')
            
        if not hasattr(self, "muc"):
            self.muc = 0.13
            logger.debug(f'muc = {self.muc}\n')
        
        if not hasattr(self, "npade"):
            self.npade = 20
            logger.debug(f'npade = {self.npade}\n')
            
        if not hasattr(self, "nstemp"):
            self.nstemp = 1
            logger.debug(f'nstemp = {self.nstemp}\n')
        else:
            self.nstemp = int(self.nstemp)
            logger.debug(f'nstemp = {self.nstemp}\n')
        
        if not hasattr(self, "temps"):
            self.temps = []
            logger.debug(f'temps = {self.temps}\n')
        else:
            self.temps = self.temps.split()
            logger.debug(f'temps = {self.temps}\n') 
            
    @classmethod
    def init_from_config(cls, config: dict):

        work_path         = config['work_path']            ; del config['work_path']
        submit_job_system = config['submit_job_system']    ; del config['submit_job_system']
        input_file_path   = config['input_file_path']      ; del config['input_file_path']
        
        self = cls(
            work_path=work_path,
            submit_job_system=submit_job_system,
            input_file_path=input_file_path,
            **config,
        )
        return self
    
    def get_hspp_for_EPW(self):
        """
        This method is to get high symmetry paths and points
        """ 
        lat     = self.ase_type.cell.get_bravais_lattice()
        pstring = lat.special_path   # GHNGPH,PN

        # 获得高对称点路径 path_name_list
        _plist  = [[ p for p in pp if not p.isdigit()] for pp in pstring.split(",")] # [['G', 'H', 'N', 'G', 'P', 'H'], ['P', 'N']]

        high_symmetry_type = np.argmax([len(_plist) for plist in _plist]) # high_symmetry_type = 0 
        path_name_list = _plist[high_symmetry_type] # ['G', 'H', 'N', 'G', 'P', 'H']
        # path_name_list = list(chain.from_iterable(_plist)) # ['G', 'H', 'N', 'G', 'P', 'H', 'P', 'N']
        print(f"the high symmetry points path: \n{path_name_list}")

        # 获得高对称路径相应的坐标
        special_points   = lat.get_special_points() # {'G': array([0, 0, 0]), 'H': array([ 0.5, -0.5,  0.5]), 'P': array([0.25, 0.25, 0.25]), 'N': array([0. , 0. , 0.5])}
        path_coords      = [list(special_points[point_name]) for point_name in path_name_list] # [[0,     0,    0], 
                                                                                               #  [0.5,  -0.5,  0.5],
                                                                                               #  [0.0,   0.0,  0.5],
                                                                                               #  [0,     0,    0], 
                                                                                               #  [0.25,  0.25, 0.25], 
                                                                                               #  [0.5,  -0.5,  0.5]]
        path_name_coords = list(zip(path_name_list, path_coords)) # [('G', [0,     0,    0]), 
                                                                  #  ('H', [0.5,  -0.5,  0.5]), 
                                                                  #  ('N', [0.0,   0.0,  0.5]), 
                                                                  #  ('G', [0,     0,    0]), 
                                                                  #  ('P', [0.25,  0.25, 0.25]), 
                                                                  #  ('H', [0.5,  -0.5,  0.5])]

        # 输出分数坐标的高对称路径
        print("Print Fractional Coordinates of Reciprocal Lattice ! ")
        for name, dirt in path_name_coords:
            print("{:<10.6f} {:<10.6f} {:<10.6f} {:<4}".format(dirt[0], dirt[1], dirt[2], name))

        # 输出倒格子晶格
        print("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
        for vector in self.reciprocal_plattice:
            print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))

        # 输出高对称路径的二维路径投影 
        print("Print projected high symmetry path")
        print("倒格子的单位是 2pi/alat")
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
        print(string_names)
        print(string_coord)

        # 输出分数坐标的高对称路径以EPW或者wannier90可以识别的方式
        path_name_coords_for_EPW = []
        for pathname_begin, pathname_end in zip(path_name_coords[:-1], path_name_coords[1:]):
            hsppinfo = f"{pathname_begin[0]:<2} {pathname_begin[1][0]:<+5.4f} {pathname_begin[1][1]:<+5.4f} {pathname_begin[1][2]:<+5.4f} {pathname_end[0]:<2} {pathname_end[1][0]:<+5.4f} {pathname_end[1][1]:<+5.4f} {pathname_end[1][2]:<+5.4f}"
            print(hsppinfo) # hsppinfo = G  0.0000 0.0000 0.0000 H  0.5000 -0.5000 0.5000
            path_name_coords_for_EPW.append(hsppinfo)

        return path_name_coords_for_EPW 

    def check_kgrid_qgird(self, grid:str):
        """
        This method is to check kgrid and qgrid
        """
        grid = grid.split()
        if len(grid) != 3:
            logger.error(f"nk must be 3 numbers, such as nk='4 4 4'")
            sys.exit(1)
        else:
            grid = [int(nk) for nk in grid]
            return grid
