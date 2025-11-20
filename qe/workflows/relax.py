"""
结构优化工作流

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from workflows.base import BaseWorkflow
from config.relax import RelaxConfig


class RelaxWorkflow(BaseWorkflow):
    """结构优化工作流"""

    def __init__(self, config: RelaxConfig):
        super().__init__(config)

    def prepare_input(self) -> Path:
        """准备relax.in输入文件"""
        input_file = self.work_dir / 'relax.in'

        with open(input_file, 'w') as f:
            # &CONTROL namelist
            f.write("&CONTROL\n")
            f.write(f"  calculation = '{self.config.get_param('calculation', 'relax')}'\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write(f"  pseudo_dir = './'\n")
            f.write(f"  outdir = './tmp/'\n")
            f.write(f"  nstep = {self.config.get_param('nstep', 100)}\n")
            f.write("/\n\n")

            # &SYSTEM namelist
            f.write("&SYSTEM\n")
            f.write(f"  ecutwfc = {self.config.get_param('ecutwfc')}\n")
            f.write(f"  ecutrho = {self.config.get_param('ecutrho')}\n")
            f.write("/\n\n")

            # &ELECTRONS namelist
            f.write("&ELECTRONS\n")
            f.write(f"  conv_thr = {self.config.get_param('conv_thr')}\n")
            f.write("/\n\n")

            # &IONS namelist (仅优化需要)
            if 'relax' in self.config.get_param('calculation'):
                f.write("&IONS\n")
                f.write("/\n\n")

            # TODO: 添加ATOMIC_SPECIES, CELL_PARAMETERS, ATOMIC_POSITIONS, K_POINTS

        return input_file

    def generate_submit_script(self) -> Path:
        """生成提交脚本"""
        script_file = self.work_dir / 's1_relax.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 结构优化计算\n\n")
            f.write("mpirun -np 8 pw.x < relax.in > relax.out\n")

        script_file.chmod(0o755)
        return script_file
