"""
VASP后处理工作流

从vasp_run.py重构而来
提供声子、电子数据后处理和文件清理功能

作者：Claude (重构自原始代码)
创建时间：2025-11-20
"""

import os
import re
import sys
import shutil
import logging
from argparse import ArgumentParser
from pathlib import Path
from itertools import chain

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

from vasp.config import config
from vasp.vasp_base import vasp_base

logger = logging.getLogger(__name__)


class PostprocessWorkflow(vasp_base):
    """
    VASP数据后处理工作流

    支持多种后处理模式：
    - dispband/dfptband: 声子能带
    - dispphdos/dfptphdos: 声子态密度
    - eband: 电子能带
    - eledos: 电子态密度
    - hspp: 高对称路径
    """

    def __init__(self, args: ArgumentParser) -> None:
        # 读取配置
        self._config = config(args).read_config()
        self.input_file_path = Path(self._config['input_file_path'])
        if self._config['work_path']:
            self.work_path = Path(self._config['work_path'])
        else:
            self.work_path = Path.cwd()
        self.mode = self._config['mode']
        self.ase_type = read(self.input_file_path)
        self.struct_type = AseAtomsAdaptor.get_structure(self.ase_type)
        self.get_struct_info(self.struct_type, self.work_path)

        if self.mode == "dispband" or self.mode == "dfptband":
            self.post_progress_phono_band()
        if self.mode == "dispphdos" or self.mode == "dfptphdos":
            self.post_progress_phono_dos()
        if self.mode == "eband":
            self.post_progress_eletron_band()
        if self.mode == "eledos":
            self.post_progress_eletron_dos()
        if self.mode == "hspp":
            self.read_hspp(self.work_path.joinpath("KPOINTS"))

    # 绘制 phonoband
    def post_progress_phono_band(self):

        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            logger.error("    you have to specify the supercell=[?,?,?]. If you didn't specify it, maybe somthing wrong will occur !")
            sys.exit(1)

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [20, 20, 20]
            logger.debug("you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if self.mode == "dispband":
            logger.debug("Run disp-band-post-progress-module")
            _disp_num = len(list(Path(self.work_path).glob("disp-*")))
            disp_num = str(_disp_num).rjust(3, '0')
            cwd = os.getcwd()
            os.chdir(self.work_path)
            logger.debug("Run the order `phonopy -f disp-{001..%s}/vasprun.xml`" % (disp_num))
            logger.debug("Please confirm  disp-{001..%s} are all correctively computed !!!" % (disp_num))
            os.system("phonopy -f disp-{001..%s}/vasprun.xml" % (disp_num))
            os.chdir(cwd)

            vaspkitflag = input("If you have installed vaspkit and you want to use it, input: Yes\n")
            if vaspkitflag:
                cwd = os.getcwd()
                os.chdir(self.work_path)
                os.system('echo -e "3\n305\n3" | vaspkit')
                shutil.copy("KPATH.phonopy", "band.conf")
                diminfo = "DIM={}".format(' '.join(list(map(str, self.supercell))))
                os.system("sed -i '2s/.*/{}/' band.conf".format(diminfo))
                os.system("sed -i '/FORCE_CONSTANTS = READ/d' band.conf")
                mpinfo = "MP={}".format(' '.join(list(map(str, self.mp))))
                os.system("sed -i '6s/.*/{}/' band.conf".format(mpinfo))
                os.chdir(cwd)
            else:
                path_name_list, path_coords = self.get_hspp(self.ase_type, self.autoselect)
                self.write_disp_band_conf(
                    self.work_path,
                    self.species,
                    self.supercell,
                    path_name_list,
                    path_coords
                )

            cwd = os.getcwd()
            os.chdir(self.work_path)
            # traditional method
            logger.debug("Run the order `phonopy -p -s band.conf -c POSCAR-init`")
            logger.debug("Please confirm the POSCAR-init exist !!!")
            os.system("phonopy -p -s band.conf -c POSCAR-init")
            logger.debug("Run the order `phonopy-bandplot  --gnuplot> band.dat`")
            logger.debug("band.dat is the data for plot phononband in Origin")
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)

        elif self.mode == "dfptband":
            logger.debug("Run dfpt-band-post-progress-module")
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy --fc vasprun.xml")
            os.chdir(cwd)

            vaspkitflag = input("If you have installed vaspkit and you want to use it, input: Yes\n")
            if vaspkitflag:
                cwd = os.getcwd()
                os.chdir(self.work_path)
                os.system('echo -e "3\n305" | vaspkit')
                shutil.copy("KPATH.phonopy", "band.conf")
                diminfo = "DIM={}".format(' '.join(list(map(str, self.supercell))))
                os.system("sed -i '2s/.*/{}/' band.conf".format(diminfo))
                os.chdir(cwd)
            else:
                special_points, path_coords = self.get_hspp(self.ase_type, self.autoselect)
                self.write_dfpt_band_conf(
                    self.work_path,
                    self.species,
                    self.supercell,
                    special_points,
                    path_coords
                )

            cwd = os.getcwd()
            os.chdir(self.work_path)
            # traditional method
            os.system("phonopy --dim='{}' -p -s band.conf -c POSCAR-init".format(' '.join(list(map(str, self.supercell)))))
            os.system("phonopy-bandplot  --gnuplot> band.dat")
            os.chdir(cwd)

    # 绘制 phonodos
    def post_progress_phono_dos(self):

        if "phdos" in self._config:
            self.pdos = self._config['pdos']
        else:
            self.pdos = "AUTO"

        if "mp" in self._config:
            _mp = self._config['mp'].split()
            self.mp = list(map(int, _mp))
        else:
            self.mp = [20, 20, 20]
            logger.debug("you didn't specify the mp='? ? ?', the program will set default mp=[8,8,8]")

        if "supercell" in self._config:
            _supercell = self._config['supercell'].split()
            self.supercell = list(map(int, _supercell))
        else:
            logger.error("    you have to specify the supercell=[?,?,?]. If you didn't specify it, maybe somthing wrong will occur !")
            sys.exit(1)

        if "tmin" in self._config:
            self.tmin = self._config['tmin']
        else:
            self.tmin = 0

        if "tmax" in self._config:
            self.tmax = self._config['tmax']
        else:
            self.tmax = 5000

        if "tstep" in self._config:
            self.tstep = self._config['tstep']
        else:
            self.tstep = 100

        if self.mode == "dispphdos":
            # 获得total_dos.dat
            self.write_disp_mesh_conf(
                self.work_path,
                self.species,
                self.supercell,
                self.mp,
                self.tmin,
                self.tmax,
                self.tstep,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf")  # -p: dos plot   -t: thermal eletron print
            os.system("phonopy -p mesh.conf -c {}".format(self.input_file_path.name))  # 获得 total_dos.dat
            os.chdir(cwd)

            # 获得pdos.dat
            self.write_disp_phdos_conf(
                self.work_path,
                self.species,
                self.supercell,
                self.mp,
                self.pdos,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p pdos.conf -c {}".format(self.input_file_path.name))  # 获得 total_dos.dat
            os.chdir(cwd)

        elif self.mode == "dfptphdos":
            # 获得total_dos.dat
            self.write_dfpt_mesh_conf(
                self.work_path,
                self.species,
                self.supercell,
                self.mp,
                self.tmin,
                self.tmax,
                self.tstep,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p -t mesh.conf")  # -p: dos plot   -t: thermal eletron print
            os.system("phonopy -p mesh.conf -c {}".format(self.input_file_path.name))  # 获得 total_dos.dat
            os.chdir(cwd)

            # 获得pdos.dat
            self.write_dfpt_phdos_conf(
                self.work_path,
                self.species,
                self.supercell,
                self.mp,
                self.pdos,
            )
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system("phonopy -p pdos.conf -c {}".format(self.input_file_path.name))  # 获得 total_dos.dat
            os.chdir(cwd)

    # 绘制 eband
    def post_progress_eletron_band(self):

        import matplotlib.pyplot as plt
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import BSPlotter

        # 检查费米能级
        vasprunxml_path = self.work_path.joinpath("vasprun.xml")
        self.check_efermi_energy(vasprunxml_path)
        vasprun = Vasprun(vasprunxml_path)
        e_fermi_fromband = vasprun.efermi
        logger.debug("Check the E-fermi( {} ) is equal to e_fermi_fromscf whether or not at last time !".format(e_fermi_fromband))

        eband = vasprun.get_band_structure(line_mode=True)

        e_fermi_fromvasp = vasprun.efermi

        # set figure parameters, draw figure
        bsplotter = BSPlotter(bs=eband)
        bsplotter.bs_plot_data()
        bsplotter.get_plot()
        ebandpng_path = self.work_path.joinpath('eband.png')
        bsplotter.save_plot(
            ebandpng_path,
            img_format='png',
            ylim=[-5, 5],
        )

    # 绘制 eledos
    def post_progress_eletron_dos(self):

        import matplotlib.pyplot as plt
        from pymatgen.io.vasp.outputs import Vasprun
        from pymatgen.electronic_structure.plotter import DosPlotter

        # 检查费米能级
        vasprunxml_path = self.work_path.joinpath("vasprun.xml")
        self.check_efermi_energy(vasprunxml_path)
        vasprun = Vasprun(vasprunxml_path)
        e_fermi_fromdos = vasprun.efermi
        logger.debug("Check the E-fermi( {} ) is equal to e_fermi_fromscf whether or not at last time !".format(e_fermi_fromdos))

        # 获得dos的数据
        try:
            eledos = vasprun.complete_dos_normalized
            logger.debug("Success getting the complete_dos_normalized")
        except:
            eledos = vasprun.complete_dos
            logger.debug("Success getting the complete_dos")

        # 处理dos的数据，并决定
        dosplotter = DosPlotter()
        self.pdostype = self._config.get('pdostype', None)
        if self.pdostype == "ele":
            dosplotter.add_dos_dict(eledos.get_element_dos())
            logger.debug("Success getting the dos projected to elements")
        elif self.pdostype == "spd":
            dosplotter.add_dos_dict(eledos.get_spd_dos())
            logger.debug("Success getting the dos projected to spd-orbits")
        elif self.pdostype == "elespd":
            dosplotter.add_dos_dict(eledos.get_element_spd_dos())
            logger.debug("Success getting the dos projected to elements and spd-orbits")
        else:
            logger.debug("NOTES: You set nothing for pdostype")
            dosplotter.add_dos_dict(eledos.get_element_dos())
            logger.debug("Default value is projected to element. Success getting the dos projected to elements")
        dosplotter.get_plot()
        eledospng_path = self.work_path.joinpath('eledos.png')
        dosplotter.save_plot(
            filename=eledospng_path,
            img_format='png',
            xlim=[-5, 5],
            ylim=[0, 20],
        )

    # 创建band.conf  目的为了获得 band of phonon
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
            path_coords = list(chain.from_iterable(path_coords))
            path_coords = list(map(str, path_coords))
            f.write("BAND={}                 \n".format(' '.join(path_coords)))

    # 创建band.conf  目的为了获得 band of phonon
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
            path_coords = list(chain.from_iterable(path_coords))
            path_coords = list(map(str, path_coords))
            f.write("BAND={}                 \n".format(' '.join(path_coords)))

    # 创建mesh.conf  目的为了获得 thermal_properties.yaml
    def write_disp_mesh_conf(
            self,
            mesh_conf_dirpath,
            species,
            supercell,
            mp,
            tmin,
            tmax,
            tstep,
    ):

        __species = [spe.name for spe in species]
        __supercell = list(map(str, supercell))
        __mp = list(map(str, mp))
        mesh_conf_filepath = os.path.join(mesh_conf_dirpath, "mesh.conf")
        with open(mesh_conf_filepath, "w") as f:
            f.write("ATOM_NAME={}            \n".format(' '.join(__species)))
            f.write("DIM={}                  \n".format(' '.join(__supercell)))
            f.write("MP ={}                  \n".format(' '.join(__mp)))
            f.write("TPROP=T                 \n")
            f.write("TMIN={}                 \n".format(tmin))
            f.write("TMAX={}                 \n".format(tmax))
            f.write("TSTEP={}                \n".format(tstep))

    # 创建mesh.conf  目的为了获得 thermal_properties.yaml
    def write_dfpt_mesh_conf(
            self,
            mesh_conf_dirpath,
            species,
            supercell,
            mp,
            tmin,
            tmax,
            tstep,
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
            f.write("TPROP=T                 \n")
            f.write("TMIN={}                 \n".format(tmin))
            f.write("TMAX={}                 \n".format(tmax))
            f.write("TSTEP={}                \n".format(tstep))

    # 创建pdos.conf  目的为了获得 pdos
    def write_disp_phdos_conf(
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
    def write_dfpt_phdos_conf(
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

    # 检查费米能级
    def check_efermi_energy(self, vasprunxml_path):

        scfoutcar_path = self.work_path.absolute().parent.joinpath("scf", "OUTCAR")
        # 检查费米能级
        e_fermi_fromscf = os.popen(f"grep E-fermi {scfoutcar_path} | tail -n 1 " + "| awk '{print $3}' ").read().strip("\n")
        e_fermi_dos_band = os.popen(f"grep efermi  {vasprunxml_path}" + "| awk '{print $3}' ").read().strip("\n")
        logger.debug("You have to confirm that the Fermi energy is from scf/OUTCAR. Because the Fermi energy in dos/DOSCAR is not accurate")
        logger.debug("You can use 'grep E-fermi scf/OUTCAR' to check the Fermi energy by yourself !")
        logger.debug("E-fermi in scf is {}".format(e_fermi_fromscf))
        logger.debug("E-fermi in dos is {}".format(e_fermi_dos_band))
        logger.debug("The program will use `e_fermi_fromscf` to cover the `e_fermi_dos_band`")
        if abs(float(e_fermi_fromscf) - float(e_fermi_dos_band)) > 0.0001:
            replace_efermi_in_vasprunxml = """ sed -E -i.bak """ + \
                                           """ 's/<i name="efermi">\s*[0-9]+\.[0-9]+\s*<\/i>/<i name="efermi">    {} <\/i>/' """.format(
                                               e_fermi_fromscf) + \
                                           """ {} """.format(vasprunxml_path)
            cwd = os.getcwd()
            os.chdir(self.work_path)
            os.system(replace_efermi_in_vasprunxml)
            os.chdir(cwd)
        logger.debug("logger.debugIf you wanna plot band or dos by yourself, you'd better replace efermi in DOSCAR with that in scf/OUTCAR")


class ClearWorkflow:
    """
    VASP文件清理工作流

    提供两种清理模式：
    - clearall: 清理所有非必要文件
    - clearopt: 清理优化过程的中间文件
    """

    def __init__(self, args: ArgumentParser) -> None:
        # 读取配置
        self._config = config(args).read_config()
        self.work_path = Path(self._config['work_path'])
        self.mode = self._config['mode']

        if self.mode == "clearall":
            self.clear_all()
        elif self.mode == "clearopt":
            self.clear_for_opt()

    def clear_all(self):
        """
        删除所有文件，除了:
        - POSCAR, PPOSCAR, POTCAR, OUTCAR
        - INCAR*
        - *.sh, *.vasp, *.slurm
        """
        reserved_files = ["POTCAR", "OUTCAR", "CONTCAR"]
        current_files = os.listdir(self.work_path)
        for file in current_files:
            if file in reserved_files:
                pass
            elif "INCAR" in file:
                pass
            elif "POSCAR" in file:
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".sh":
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".vasp":
                pass
            elif Path(self.work_path).joinpath(file).suffix == ".slurm":
                pass
            else:
                os.system(f"rm -fr {str(Path(self.work_path).joinpath(file))}")

    def clear_for_opt(self):
        """清理优化过程的中间文件"""
        files = ["CHG", "CHGCAR", "DOSCAR", "EIGENVAL", "FERMI_ENERGY",
                 "OSZICAR", "PCDAT", "PROCAR", "REPORT", "WAVECAR",
                 "XDATCAR", "vasprun.xml"]
        for file in files:
            os.system(f"rm -fr {str(Path(self.work_path).joinpath(file))}")
