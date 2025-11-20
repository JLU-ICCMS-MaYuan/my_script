"""
自洽计算工作流

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import Optional, Dict, Any
import logging

from workflows.base import BaseWorkflow
from config.base import BaseConfig
from utils.structure import StructureReader

logger = logging.getLogger(__name__)


class ScfWorkflow(BaseWorkflow):
    """自洽计算（SCF）工作流"""

    def __init__(self, config: BaseConfig, structure_file: Optional[Path] = None):
        """
        初始化SCF工作流

        Parameters
        ----------
        config : BaseConfig
            配置对象
        structure_file : Path, optional
            结构文件路径（POSCAR/cif等）
        """
        super().__init__(config)
        self.structure_file = structure_file

        if structure_file:
            self.structure = StructureReader.read_structure(structure_file)
        else:
            self.structure = None

    def prepare_input(self) -> Path:
        """准备scf.in输入文件"""
        input_file = self.work_dir / 'scf.in'

        with open(input_file, 'w') as f:
            # &CONTROL namelist
            f.write("&CONTROL\n")
            f.write("  calculation = 'scf'\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  pseudo_dir = './'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  verbosity = 'high'\n")

            # 可选参数
            if self.config.get_param('tstress'):
                f.write("  tstress = .true.\n")
            if self.config.get_param('tprnfor'):
                f.write("  tprnfor = .true.\n")

            f.write("/\n\n")

            # &SYSTEM namelist
            f.write("&SYSTEM\n")

            # 基本参数
            if self.structure:
                f.write(f"  ibrav = 0\n")
                f.write(f"  nat = {len(self.structure['species'])}\n")

                # 计算元素种类
                unique_species = set(self.structure['species'])
                f.write(f"  ntyp = {len(unique_species)}\n")

            # 能量截断
            ecutwfc = self.config.get_param('ecutwfc')
            ecutrho = self.config.get_param('ecutrho')
            if ecutwfc:
                f.write(f"  ecutwfc = {ecutwfc}\n")
            if ecutrho:
                f.write(f"  ecutrho = {ecutrho}\n")

            # 占据方式
            occupations = self.config.get_param('occupations', 'smearing')
            f.write(f"  occupations = '{occupations}'\n")

            if occupations == 'smearing':
                smearing = self.config.get_param('smearing', 'mp')
                degauss = self.config.get_param('degauss', 0.02)
                f.write(f"  smearing = '{smearing}'\n")
                f.write(f"  degauss = {degauss}\n")

            # 自旋相关
            if self.config.get_param('nspin'):
                f.write(f"  nspin = {self.config.get_param('nspin')}\n")

            # vdW修正
            vdw_corr = self.config.get_param('vdw_corr')
            if vdw_corr:
                f.write(f"  vdw_corr = '{vdw_corr}'\n")

            f.write("/\n\n")

            # &ELECTRONS namelist
            f.write("&ELECTRONS\n")

            conv_thr = self.config.get_param('conv_thr', '1.0d-8')
            f.write(f"  conv_thr = {conv_thr}\n")

            mixing_beta = self.config.get_param('mixing_beta', 0.3)
            f.write(f"  mixing_beta = {mixing_beta}\n")

            electron_maxstep = self.config.get_param('electron_maxstep', 100)
            f.write(f"  electron_maxstep = {electron_maxstep}\n")

            f.write("/\n\n")

            # 如果有结构信息，写入原子信息
            if self.structure:
                self._write_atomic_info(f)

        logger.info(f"SCF输入文件已生成: {input_file}")
        return input_file

    def _write_atomic_info(self, f):
        """写入原子种类、晶格、原子位置和K点信息"""

        # ATOMIC_SPECIES
        f.write("ATOMIC_SPECIES\n")
        unique_species = sorted(set(self.structure['species']))

        # 获取赝势信息
        pseudopotentials = self.config.get_param('pseudopotentials', {})

        for species in unique_species:
            # 从元素周期表获取质量（简化版）
            mass = self._get_atomic_mass(species)
            pseudo = pseudopotentials.get(species, f"{species}.UPF")
            f.write(f"  {species}  {mass:.3f}  {pseudo}\n")
        f.write("\n")

        # CELL_PARAMETERS
        f.write("CELL_PARAMETERS angstrom\n")
        lattice = self.structure['lattice']
        for row in lattice:
            f.write(f"  {row[0]:15.10f} {row[1]:15.10f} {row[2]:15.10f}\n")
        f.write("\n")

        # ATOMIC_POSITIONS
        coord_type = self.config.get_param('coord_type', 'crystal')
        f.write(f"ATOMIC_POSITIONS {coord_type}\n")

        for species, pos in zip(self.structure['species'], self.structure['positions']):
            f.write(f"  {species:3s} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n")
        f.write("\n")

        # K_POINTS
        kpoints = self.config.get_param('kpoints', (8, 8, 8))
        kshift = self.config.get_param('kshift', (0, 0, 0))

        f.write("K_POINTS automatic\n")
        f.write(f"  {kpoints[0]} {kpoints[1]} {kpoints[2]} {kshift[0]} {kshift[1]} {kshift[2]}\n")
        f.write("\n")

    def _get_atomic_mass(self, element: str) -> float:
        """获取元素原子质量（简化版）"""
        # 常见元素的原子质量
        masses = {
            'H': 1.008, 'He': 4.003, 'Li': 6.941, 'C': 12.011, 'N': 14.007,
            'O': 15.999, 'F': 18.998, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982,
            'Si': 28.086, 'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'K': 39.098,
            'Ca': 40.078, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
            'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
            'Ga': 69.723, 'Ge': 72.64, 'As': 74.922, 'Se': 78.96, 'Br': 79.904,
            'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.96,
            'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710,
            'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Ba': 137.327, 'La': 138.905,
            'W': 183.84, 'Pt': 195.084, 'Au': 196.967, 'Pb': 207.2, 'Bi': 208.980,
        }
        return masses.get(element, 1.0)

    def generate_submit_script(self) -> Path:
        """生成SCF计算提交脚本"""
        script_file = self.work_dir / 's_scf.sh'

        # 获取并行核数
        nprocs = self.config.get_param('nprocs', 8)
        qe_bin = self.config.get_param('qe_bin', 'pw.x')

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 自洽计算（SCF）\n\n")
            f.write(f"mpirun -np {nprocs} {qe_bin} < scf.in > scf.out\n\n")
            f.write("# 检查是否收敛\n")
            f.write("if grep -q 'convergence NOT achieved' scf.out; then\n")
            f.write("    echo '错误: SCF未收敛！'\n")
            f.write("    exit 1\n")
            f.write("elif grep -q 'convergence has been achieved' scf.out; then\n")
            f.write("    echo 'SCF计算成功收敛'\n")
            f.write("else\n")
            f.write("    echo '警告: 无法确定SCF收敛状态'\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"SCF提交脚本已生成: {script_file}")
        return script_file
