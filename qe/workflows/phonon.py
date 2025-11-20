"""
声子计算工作流

支持多种声子计算模式：
- nosplit: 不分割q点，一次性计算
- split_dyn0: 分割q点，从q=0开始
- split_assignQ: 分割q点，指定特定q点
- merge: 合并已分割的动力学矩阵
- matdyn: 计算声子频率色散
- phonodos: 计算声子态密度
- el_ph: 电声耦合计算

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import Optional, List
import logging

from workflows.base import BaseWorkflow
from config.phonon import PhononConfig

logger = logging.getLogger(__name__)


class PhononWorkflow(BaseWorkflow):
    """声子计算工作流"""

    def __init__(self, config: PhononConfig, mode: str = 'nosplit'):
        """
        初始化声子工作流

        Parameters
        ----------
        config : PhononConfig
            声子配置对象
        mode : str
            计算模式: nosplit, split_dyn0, split_assignQ, merge, matdyn, phonodos, el_ph
        """
        super().__init__(config)
        self.mode = mode

    def prepare_input(self) -> Path:
        """根据mode准备相应的输入文件"""
        if self.mode == 'nosplit':
            return self._prepare_ph_nosplit()
        elif self.mode == 'split_dyn0':
            return self._prepare_ph_split_dyn0()
        elif self.mode == 'split_assignQ':
            return self._prepare_ph_assignQ()
        elif self.mode == 'merge':
            return self._prepare_dynmat_merge()
        elif self.mode == 'matdyn':
            return self._prepare_matdyn()
        elif self.mode == 'phonodos':
            return self._prepare_phonon_dos()
        elif self.mode == 'el_ph':
            return self._prepare_el_ph()
        else:
            raise ValueError(f"未知的声子计算模式: {self.mode}")

    def _prepare_ph_nosplit(self) -> Path:
        """准备不分割q点的ph.x输入"""
        input_file = self.work_dir / 'ph.in'

        with open(input_file, 'w') as f:
            f.write("&INPUTPH\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write(f"  fildyn = 'dyn'\n")

            # 收敛参数
            tr2_ph = self.config.get_param('tr2_ph', '1.0d-14')
            f.write(f"  tr2_ph = {tr2_ph}\n")

            # 线性响应参数
            alpha_mix = self.config.get_param('alpha_mix', 0.3)
            f.write(f"  alpha_mix(1) = {alpha_mix}\n")

            # 是否计算电声耦合
            if self.config.get_param('el_ph_calculation'):
                f.write("  ldisp = .true.\n")
                f.write("  trans = .true.\n")
                f.write("  electron_phonon = 'interpolated'\n")

                sigma = self.config.get_param('el_ph_sigma', 0.005)
                nsigma = self.config.get_param('el_ph_nsigma', 10)
                f.write(f"  el_ph_sigma = {sigma}\n")
                f.write(f"  el_ph_nsigma = {nsigma}\n")
            else:
                f.write("  ldisp = .true.\n")
                f.write("  trans = .true.\n")

            f.write("/\n")

            # Q点网格
            qpoints = self.config.get_param('qpoints', '6 6 6')
            f.write(f"{qpoints}\n")

        logger.info(f"声子输入文件已生成（nosplit模式）: {input_file}")
        return input_file

    def _prepare_ph_split_dyn0(self) -> Path:
        """准备从dyn0开始分割计算的ph.x输入"""
        input_file = self.work_dir / 'ph.in'

        with open(input_file, 'w') as f:
            f.write("&INPUTPH\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write(f"  fildyn = 'dyn'\n")

            tr2_ph = self.config.get_param('tr2_ph', '1.0d-14')
            f.write(f"  tr2_ph = {tr2_ph}\n")

            alpha_mix = self.config.get_param('alpha_mix', 0.3)
            f.write(f"  alpha_mix(1) = {alpha_mix}\n")

            # 分割参数
            f.write("  ldisp = .true.\n")
            f.write("  trans = .true.\n")
            f.write("  start_q = 1\n")
            f.write("  last_q = 1\n")  # 只计算q=0点

            f.write("/\n")

            qpoints = self.config.get_param('qpoints', '6 6 6')
            f.write(f"{qpoints}\n")

        logger.info(f"声子输入文件已生成（split_dyn0模式）: {input_file}")
        return input_file

    def _prepare_ph_assignQ(self) -> Path:
        """准备指定q点计算的ph.x输入"""
        input_file = self.work_dir / 'ph.in'

        start_q = self.config.get_param('start_q', 2)
        last_q = self.config.get_param('last_q', 2)

        with open(input_file, 'w') as f:
            f.write("&INPUTPH\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write(f"  fildyn = 'dyn'\n")

            tr2_ph = self.config.get_param('tr2_ph', '1.0d-14')
            f.write(f"  tr2_ph = {tr2_ph}\n")

            alpha_mix = self.config.get_param('alpha_mix', 0.3)
            f.write(f"  alpha_mix(1) = {alpha_mix}\n")

            f.write("  ldisp = .true.\n")
            f.write("  trans = .true.\n")
            f.write(f"  start_q = {start_q}\n")
            f.write(f"  last_q = {last_q}\n")
            f.write("  recover = .true.\n")  # 从之前的计算恢复

            f.write("/\n")

            qpoints = self.config.get_param('qpoints', '6 6 6')
            f.write(f"{qpoints}\n")

        logger.info(f"声子输入文件已生成（assignQ模式 q={start_q}-{last_q}）: {input_file}")
        return input_file

    def _prepare_dynmat_merge(self) -> Path:
        """准备动力学矩阵合并的q2r.x输入"""
        input_file = self.work_dir / 'q2r.in'

        with open(input_file, 'w') as f:
            f.write("&INPUT\n")
            f.write("  fildyn = 'dyn'\n")
            f.write("  flfrc = 'force_constants.fc'\n")
            f.write("/\n")

        logger.info(f"q2r输入文件已生成: {input_file}")
        return input_file

    def _prepare_matdyn(self) -> Path:
        """准备matdyn.x计算声子频率色散的输入"""
        input_file = self.work_dir / 'matdyn.in'

        with open(input_file, 'w') as f:
            f.write("&INPUT\n")
            f.write("  flfrc = 'force_constants.fc'\n")
            f.write("  flfrq = 'matdyn.freq'\n")
            f.write("  fleig = 'matdyn.eig'\n")
            f.write("  asr = 'crystal'\n")

            # 可选：声学和求和规则
            asr_type = self.config.get_param('asr', 'crystal')
            f.write(f"  asr = '{asr_type}'\n")

            f.write("/\n")

            # 高对称路径点数
            npts = self.config.get_param('band_npts', 100)
            f.write(f"{npts}\n")

            # 高对称点路径（这里需要从config读取或自动生成）
            kpath = self.config.get_param('kpath', [])
            if kpath:
                for point in kpath:
                    f.write(f"{point[0]:10.6f} {point[1]:10.6f} {point[2]:10.6f}  {point[3]}\n")
            else:
                # 默认Gamma-X路径
                f.write("0.0  0.0  0.0  50  ! Gamma\n")
                f.write("0.5  0.0  0.5  50  ! X\n")

        logger.info(f"matdyn输入文件已生成: {input_file}")
        return input_file

    def _prepare_phonon_dos(self) -> Path:
        """准备matdyn.x计算声子态密度的输入"""
        input_file = self.work_dir / 'phdos.in'

        with open(input_file, 'w') as f:
            f.write("&INPUT\n")
            f.write("  flfrc = 'force_constants.fc'\n")
            f.write("  flfrq = 'phonon.dos'\n")
            f.write("  dos = .true.\n")

            # DOS相关参数
            dos_emin = self.config.get_param('dos_emin', 0.0)
            dos_emax = self.config.get_param('dos_emax', 1000.0)
            dos_ndos = self.config.get_param('dos_ndos', 1000)
            dos_delta = self.config.get_param('dos_delta', 1.0)

            f.write(f"  dos_emin = {dos_emin}\n")
            f.write(f"  dos_emax = {dos_emax}\n")
            f.write(f"  dos_ndos = {dos_ndos}\n")
            f.write(f"  dos_delta = {dos_delta}\n")

            asr_type = self.config.get_param('asr', 'crystal')
            f.write(f"  asr = '{asr_type}'\n")

            f.write("/\n")

            # Q点网格（用于DOS）
            dos_qgrid = self.config.get_param('dos_qgrid', (20, 20, 20))
            f.write(f"{dos_qgrid[0]} {dos_qgrid[1]} {dos_qgrid[2]}\n")

        logger.info(f"声子DOS输入文件已生成: {input_file}")
        return input_file

    def _prepare_el_ph(self) -> Path:
        """准备电声耦合计算的ph.x输入"""
        input_file = self.work_dir / 'elph.in'

        with open(input_file, 'w') as f:
            f.write("&INPUTPH\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write(f"  fildyn = 'dyn'\n")

            tr2_ph = self.config.get_param('tr2_ph', '1.0d-14')
            f.write(f"  tr2_ph = {tr2_ph}\n")

            f.write("  ldisp = .true.\n")
            f.write("  trans = .true.\n")
            f.write("  electron_phonon = 'interpolated'\n")

            # 电声耦合参数
            sigma = self.config.get_param('el_ph_sigma', 0.005)
            nsigma = self.config.get_param('el_ph_nsigma', 10)
            f.write(f"  el_ph_sigma = {sigma}\n")
            f.write(f"  el_ph_nsigma = {nsigma}\n")

            f.write("/\n")

            qpoints = self.config.get_param('qpoints', '6 6 6')
            f.write(f"{qpoints}\n")

        logger.info(f"电声耦合输入文件已生成: {input_file}")
        return input_file

    def generate_submit_script(self) -> Path:
        """生成声子计算提交脚本"""
        nprocs = self.config.get_param('nprocs', 8)

        if self.mode == 'nosplit':
            return self._generate_ph_script(nprocs)
        elif self.mode in ['split_dyn0', 'split_assignQ']:
            return self._generate_ph_script(nprocs)
        elif self.mode == 'merge':
            return self._generate_q2r_script(nprocs)
        elif self.mode == 'matdyn':
            return self._generate_matdyn_script(nprocs)
        elif self.mode == 'phonodos':
            return self._generate_phdos_script(nprocs)
        elif self.mode == 'el_ph':
            return self._generate_elph_script(nprocs)
        else:
            raise ValueError(f"未知的声子计算模式: {self.mode}")

    def _generate_ph_script(self, nprocs: int) -> Path:
        """生成ph.x脚本"""
        script_file = self.work_dir / f's_phonon_{self.mode}.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write(f"# 声子计算 (模式: {self.mode})\n\n")
            f.write(f"mpirun -np {nprocs} ph.x < ph.in > ph.out\n\n")
            f.write("# 检查计算状态\n")
            f.write("if grep -q 'JOB DONE' ph.out; then\n")
            f.write("    echo '声子计算成功完成'\n")
            f.write("else\n")
            f.write("    echo '错误: 声子计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"声子计算脚本已生成: {script_file}")
        return script_file

    def _generate_q2r_script(self, nprocs: int) -> Path:
        """生成q2r.x脚本"""
        script_file = self.work_dir / 's_q2r.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 动力学矩阵转力常数\n\n")
            f.write("q2r.x < q2r.in > q2r.out\n\n")
            f.write("if [ -f force_constants.fc ]; then\n")
            f.write("    echo '力常数文件已生成'\n")
            f.write("else\n")
            f.write("    echo '错误: 力常数文件生成失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"q2r脚本已生成: {script_file}")
        return script_file

    def _generate_matdyn_script(self, nprocs: int) -> Path:
        """生成matdyn.x脚本"""
        script_file = self.work_dir / 's_matdyn.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 计算声子频率色散\n\n")
            f.write("matdyn.x < matdyn.in > matdyn.out\n\n")
            f.write("if [ -f matdyn.freq ]; then\n")
            f.write("    echo '声子频率已计算'\n")
            f.write("else\n")
            f.write("    echo '错误: 声子频率计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"matdyn脚本已生成: {script_file}")
        return script_file

    def _generate_phdos_script(self, nprocs: int) -> Path:
        """生成声子DOS计算脚本"""
        script_file = self.work_dir / 's_phdos.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 计算声子态密度\n\n")
            f.write("matdyn.x < phdos.in > phdos.out\n\n")
            f.write("if [ -f phonon.dos ]; then\n")
            f.write("    echo '声子DOS已计算'\n")
            f.write("else\n")
            f.write("    echo '错误: 声子DOS计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"声子DOS脚本已生成: {script_file}")
        return script_file

    def _generate_elph_script(self, nprocs: int) -> Path:
        """生成电声耦合计算脚本"""
        script_file = self.work_dir / 's_elph.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 电声耦合计算\n\n")
            f.write(f"mpirun -np {nprocs} ph.x < elph.in > elph.out\n\n")
            f.write("# 检查a2F文件是否生成\n")
            f.write("if [ -f _ph0/*.a2Fsave ]; then\n")
            f.write("    echo '电声耦合参数已计算'\n")
            f.write("else\n")
            f.write("    echo '警告: 未找到a2F文件'\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"电声耦合脚本已生成: {script_file}")
        return script_file
