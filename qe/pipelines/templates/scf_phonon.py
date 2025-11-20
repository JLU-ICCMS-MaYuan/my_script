"""
SCF + 声子计算 Pipeline

工作流：
1. SCF - 自洽计算
2. Phonon (nosplit) - 声子计算（不分割q点）
3. Q2R - 动力学矩阵→力常数
4. Matdyn - 声子频率色散
5. PhononDOS - 声子态密度

适用场景：
- 已优化结构的声子谱计算
- 快速声子性质评估
- 虚频检查

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import List, Dict, Any
import logging

from pipelines.base import BasePipeline

logger = logging.getLogger(__name__)


class ScfPhononPipeline(BasePipeline):
    """SCF + 声子计算Pipeline"""

    def define_steps(self) -> List[Dict[str, Any]]:
        """
        定义计算步骤

        Returns
        -------
        List[Dict[str, Any]]
            步骤定义列表
        """
        return [
            {
                'name': 'scf',
                'workflow': 'ScfWorkflow',
                'description': 'SCF自洽计算',
                'params': {
                    'occupations': 'smearing',
                    'smearing': 'mp',
                    'degauss': 0.02,
                    'tstress': True,
                    'tprnfor': True,
                },
            },
            {
                'name': 'phonon',
                'workflow': 'PhononWorkflow',
                'mode': 'nosplit',
                'description': '声子计算（不分割）',
                'params': {
                    'tr2_ph': '1.0d-14',
                    'alpha_mix': 0.3,
                },
                'depends_on': ['scf'],
            },
            {
                'name': 'q2r',
                'workflow': 'PhononWorkflow',
                'mode': 'merge',
                'description': '动力学矩阵→力常数',
                'depends_on': ['phonon'],
            },
            {
                'name': 'matdyn',
                'workflow': 'PhononWorkflow',
                'mode': 'matdyn',
                'description': '声子频率色散',
                'depends_on': ['q2r'],
            },
            {
                'name': 'phonon_dos',
                'workflow': 'PhononWorkflow',
                'mode': 'phonodos',
                'description': '声子态密度',
                'params': {
                    'dos_qgrid': (20, 20, 20),
                    'dos_emin': 0.0,
                    'dos_emax': 1000.0,
                    'dos_ndos': 1000,
                },
                'depends_on': ['q2r'],
            },
        ]

    def get_default_config(self) -> Dict[str, Any]:
        """
        获取默认配置

        Returns
        -------
        Dict[str, Any]
            默认配置参数
        """
        return {
            # 基本参数
            'ecutwfc': 80,
            'ecutrho': 640,
            'kpoints': (8, 8, 8),
            'qpoints': '6 6 6',

            # 收敛参数
            'conv_thr': '1.0d-8',
            'mixing_beta': 0.3,

            # 赝势
            'pseudopotentials': {},

            # 声子路径（立方系统）
            'kpath': [
                [0.0, 0.0, 0.0, 50, 'Γ'],
                [0.5, 0.0, 0.5, 50, 'X'],
                [0.5, 0.5, 0.5, 50, 'M'],
                [0.0, 0.0, 0.0, 1, 'Γ'],
            ],

            # 并行设置
            'nprocs': 8,
        }


class ScfPhononElphPipeline(ScfPhononPipeline):
    """SCF + 声子 + 电声耦合 Pipeline"""

    def define_steps(self) -> List[Dict[str, Any]]:
        """
        扩展：添加电声耦合计算

        Returns
        -------
        List[Dict[str, Any]]
            步骤定义列表
        """
        steps = super().define_steps()

        # 在声子计算后添加电声耦合
        elph_step = {
            'name': 'el_ph',
            'workflow': 'PhononWorkflow',
            'mode': 'el_ph',
            'description': '电声耦合计算',
            'params': {
                'el_ph_sigma': 0.005,
                'el_ph_nsigma': 10,
            },
            'depends_on': ['scf'],
        }

        # 插入到phonon之后
        steps.insert(2, elph_step)

        return steps


def example_usage():
    """使用示例"""
    from pathlib import Path

    # 已优化的结构文件
    structures = [
        Path("/path/to/optimized_structure.vasp"),
    ]

    # 配置
    config = {
        'ecutwfc': 80,
        'ecutrho': 640,
        'kpoints': (8, 8, 8),
        'qpoints': '6 6 6',
        'pseudopotentials': {
            'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
            'S': 'S.pbe-n-rrkjus_psl.1.0.0.UPF',
        },
        'nprocs': 8,
    }

    # 创建Pipeline（带电声耦合）
    pipeline = ScfPhononElphPipeline(
        structures=structures,
        work_dir=Path("./phonon_calculations"),
        config=config,
    )

    # 运行
    pipeline.run(max_workers=2)


if __name__ == "__main__":
    example_usage()
