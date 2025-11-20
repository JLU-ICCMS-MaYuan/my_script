"""
电子性质计算工作流

支持多种电子性质计算：
- bands: 能带结构计算
- dos: 态密度计算
- pdos: 投影态密度计算
- combined: 能带+DOS组合计算

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import Optional, List, Tuple
import logging

from workflows.base import BaseWorkflow
from config.electron import ElectronConfig

logger = logging.getLogger(__name__)


class ElectronWorkflow(BaseWorkflow):
    """电子性质计算工作流"""

    def __init__(self, config: ElectronConfig, mode: str = 'bands'):
        """
        初始化电子工作流

        Parameters
        ----------
        config : ElectronConfig
            电子性质配置对象
        mode : str
            计算模式: bands, dos, pdos, combined
        """
        super().__init__(config)
        self.mode = mode

    def prepare_input(self) -> Path:
        """根据mode准备相应的输入文件"""
        if self.mode == 'bands':
            return self._prepare_bands()
        elif self.mode == 'dos':
            return self._prepare_dos()
        elif self.mode == 'pdos':
            return self._prepare_pdos()
        elif self.mode == 'combined':
            # 组合模式会生成多个文件
            self._prepare_bands()
            return self._prepare_dos()
        else:
            raise ValueError(f"未知的电子计算模式: {self.mode}")

    def _prepare_bands(self) -> Path:
        """准备能带计算输入 (bands.x)"""
        # 首先需要pw.x的nscf计算
        nscf_file = self._prepare_nscf_for_bands()

        # 然后准备bands.x的输入
        bands_file = self.work_dir / 'bands.in'

        with open(bands_file, 'w') as f:
            f.write("&BANDS\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  filband = 'bands.dat'\n")

            # 自旋极化
            if self.config.get_param('lsym'):
                f.write("  lsym = .true.\n")

            f.write("/\n")

        logger.info(f"能带计算输入已生成: {bands_file}")
        return bands_file

    def _prepare_nscf_for_bands(self) -> Path:
        """准备能带计算的nscf输入"""
        nscf_file = self.work_dir / 'nscf_bands.in'

        with open(nscf_file, 'w') as f:
            f.write("&CONTROL\n")
            f.write("  calculation = 'bands'\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  pseudo_dir = './'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  verbosity = 'high'\n")
            f.write("/\n\n")

            f.write("&SYSTEM\n")
            # 这些参数应该从scf.in继承
            f.write("  ibrav = 0\n")

            # 能量截断（从config读取）
            ecutwfc = self.config.get_param('ecutwfc')
            ecutrho = self.config.get_param('ecutrho')
            if ecutwfc:
                f.write(f"  ecutwfc = {ecutwfc}\n")
            if ecutrho:
                f.write(f"  ecutrho = {ecutrho}\n")

            # 能带数（应多于占据能带）
            nbnd = self.config.get_param('nbnd')
            if nbnd:
                f.write(f"  nbnd = {nbnd}\n")

            f.write("  occupations = 'tetrahedra'\n")
            f.write("/\n\n")

            f.write("&ELECTRONS\n")
            f.write("  conv_thr = 1.0d-10\n")
            f.write("/\n\n")

            # K点路径（高对称点）
            f.write("K_POINTS crystal_b\n")
            kpath = self.config.get_param('kpath', [])

            if kpath:
                f.write(f"{len(kpath)}\n")
                for point in kpath:
                    # 格式: kx ky kz npoints
                    f.write(f"{point[0]:10.6f} {point[1]:10.6f} {point[2]:10.6f} {point[3]}\n")
            else:
                # 默认Gamma-X-M-Gamma路径（立方晶系）
                f.write("4\n")
                f.write("0.0  0.0  0.0  50  ! Gamma\n")
                f.write("0.5  0.0  0.5  50  ! X\n")
                f.write("0.5  0.5  0.5  50  ! M\n")
                f.write("0.0  0.0  0.0  1   ! Gamma\n")

        logger.info(f"能带nscf输入已生成: {nscf_file}")
        return nscf_file

    def _prepare_dos(self) -> Path:
        """准备态密度计算输入"""
        # 首先需要pw.x的nscf计算（密集k点）
        nscf_file = self._prepare_nscf_for_dos()

        # 然后准备dos.x的输入
        dos_file = self.work_dir / 'dos.in'

        with open(dos_file, 'w') as f:
            f.write("&DOS\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  fildos = 'dos.dat'\n")

            # DOS能量范围
            emin = self.config.get_param('dos_emin', -20.0)
            emax = self.config.get_param('dos_emax', 20.0)
            deltae = self.config.get_param('dos_deltae', 0.01)

            f.write(f"  Emin = {emin}\n")
            f.write(f"  Emax = {emax}\n")
            f.write(f"  DeltaE = {deltae}\n")

            f.write("/\n")

        logger.info(f"DOS计算输入已生成: {dos_file}")
        return dos_file

    def _prepare_nscf_for_dos(self) -> Path:
        """准备DOS计算的nscf输入"""
        nscf_file = self.work_dir / 'nscf_dos.in'

        with open(nscf_file, 'w') as f:
            f.write("&CONTROL\n")
            f.write("  calculation = 'nscf'\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  pseudo_dir = './'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("/\n\n")

            f.write("&SYSTEM\n")
            f.write("  ibrav = 0\n")

            ecutwfc = self.config.get_param('ecutwfc')
            ecutrho = self.config.get_param('ecutrho')
            if ecutwfc:
                f.write(f"  ecutwfc = {ecutwfc}\n")
            if ecutrho:
                f.write(f"  ecutrho = {ecutrho}\n")

            # DOS需要更多能带
            nbnd = self.config.get_param('nbnd')
            if nbnd:
                f.write(f"  nbnd = {nbnd}\n")

            f.write("  occupations = 'tetrahedra'\n")
            f.write("/\n\n")

            f.write("&ELECTRONS\n")
            f.write("  conv_thr = 1.0d-10\n")
            f.write("/\n\n")

            # DOS需要密集的k点
            f.write("K_POINTS automatic\n")
            kpoints_dense = self.config.get_param('kpoints_dense', (16, 16, 16))
            f.write(f"{kpoints_dense[0]} {kpoints_dense[1]} {kpoints_dense[2]} 0 0 0\n")

        logger.info(f"DOS nscf输入已生成: {nscf_file}")
        return nscf_file

    def _prepare_pdos(self) -> Path:
        """准备投影态密度计算输入 (projwfc.x)"""
        pdos_file = self.work_dir / 'pdos.in'

        with open(pdos_file, 'w') as f:
            f.write("&PROJWFC\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  filpdos = 'pdos'\n")

            # 能量范围
            emin = self.config.get_param('dos_emin', -20.0)
            emax = self.config.get_param('dos_emax', 20.0)
            deltae = self.config.get_param('dos_deltae', 0.01)

            f.write(f"  Emin = {emin}\n")
            f.write(f"  Emax = {emax}\n")
            f.write(f"  DeltaE = {deltae}\n")

            # 展宽
            degauss = self.config.get_param('pdos_degauss', 0.01)
            f.write(f"  degauss = {degauss}\n")

            f.write("/\n")

        logger.info(f"PDOS计算输入已生成: {pdos_file}")
        return pdos_file

    def generate_submit_script(self) -> Path:
        """生成电子性质计算提交脚本"""
        nprocs = self.config.get_param('nprocs', 8)

        if self.mode == 'bands':
            return self._generate_bands_script(nprocs)
        elif self.mode == 'dos':
            return self._generate_dos_script(nprocs)
        elif self.mode == 'pdos':
            return self._generate_pdos_script(nprocs)
        elif self.mode == 'combined':
            return self._generate_combined_script(nprocs)
        else:
            raise ValueError(f"未知的电子计算模式: {self.mode}")

    def _generate_bands_script(self, nprocs: int) -> Path:
        """生成能带计算脚本"""
        script_file = self.work_dir / 's_bands.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 能带结构计算\n\n")

            f.write("echo '步骤1: NSCF计算（能带路径）'\n")
            f.write(f"mpirun -np {nprocs} pw.x < nscf_bands.in > nscf_bands.out\n\n")

            f.write("if ! grep -q 'convergence has been achieved' nscf_bands.out; then\n")
            f.write("    echo '错误: NSCF计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("echo '步骤2: 提取能带'\n")
            f.write("bands.x < bands.in > bands.out\n\n")

            f.write("if [ -f bands.dat ]; then\n")
            f.write("    echo '能带数据已生成'\n")
            f.write("else\n")
            f.write("    echo '错误: 能带数据生成失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"能带计算脚本已生成: {script_file}")
        return script_file

    def _generate_dos_script(self, nprocs: int) -> Path:
        """生成DOS计算脚本"""
        script_file = self.work_dir / 's_dos.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 态密度计算\n\n")

            f.write("echo '步骤1: NSCF计算（密集k点）'\n")
            f.write(f"mpirun -np {nprocs} pw.x < nscf_dos.in > nscf_dos.out\n\n")

            f.write("if ! grep -q 'convergence has been achieved' nscf_dos.out; then\n")
            f.write("    echo '错误: NSCF计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("echo '步骤2: 计算DOS'\n")
            f.write("dos.x < dos.in > dos.out\n\n")

            f.write("if [ -f dos.dat ]; then\n")
            f.write("    echo 'DOS数据已生成'\n")
            f.write("else\n")
            f.write("    echo '错误: DOS数据生成失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"DOS计算脚本已生成: {script_file}")
        return script_file

    def _generate_pdos_script(self, nprocs: int) -> Path:
        """生成PDOS计算脚本"""
        script_file = self.work_dir / 's_pdos.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 投影态密度计算\n\n")

            f.write("echo '步骤1: NSCF计算（密集k点）'\n")
            f.write(f"mpirun -np {nprocs} pw.x < nscf_dos.in > nscf_dos.out\n\n")

            f.write("if ! grep -q 'convergence has been achieved' nscf_dos.out; then\n")
            f.write("    echo '错误: NSCF计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("echo '步骤2: 计算PDOS'\n")
            f.write(f"mpirun -np {nprocs} projwfc.x < pdos.in > pdos.out\n\n")

            f.write("if ls pdos*.dat 1> /dev/null 2>&1; then\n")
            f.write("    echo 'PDOS数据已生成'\n")
            f.write("else\n")
            f.write("    echo '错误: PDOS数据生成失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"PDOS计算脚本已生成: {script_file}")
        return script_file

    def _generate_combined_script(self, nprocs: int) -> Path:
        """生成能带+DOS组合计算脚本"""
        script_file = self.work_dir / 's_bands_dos.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 能带+DOS组合计算\n\n")

            # 能带部分
            f.write("echo '========== 能带计算 =========='\n")
            f.write("echo '步骤1: NSCF计算（能带路径）'\n")
            f.write(f"mpirun -np {nprocs} pw.x < nscf_bands.in > nscf_bands.out\n\n")

            f.write("if ! grep -q 'convergence has been achieved' nscf_bands.out; then\n")
            f.write("    echo '错误: 能带NSCF计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("echo '步骤2: 提取能带'\n")
            f.write("bands.x < bands.in > bands.out\n\n")

            # DOS部分
            f.write("echo '\\n========== DOS计算 =========='\n")
            f.write("echo '步骤3: NSCF计算（密集k点）'\n")
            f.write(f"mpirun -np {nprocs} pw.x < nscf_dos.in > nscf_dos.out\n\n")

            f.write("if ! grep -q 'convergence has been achieved' nscf_dos.out; then\n")
            f.write("    echo '错误: DOS NSCF计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("echo '步骤4: 计算DOS'\n")
            f.write("dos.x < dos.in > dos.out\n\n")

            f.write("echo '\\n========== 计算完成 =========='\n")
            f.write("if [ -f bands.dat ] && [ -f dos.dat ]; then\n")
            f.write("    echo '能带和DOS数据已全部生成'\n")
            f.write("else\n")
            f.write("    echo '错误: 部分数据生成失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"能带+DOS组合脚本已生成: {script_file}")
        return script_file
