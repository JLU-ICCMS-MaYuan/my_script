"""
结构优化 + 电子性质计算 Pipeline

工作流：
1. Relax - 结构优化
2. SCF - 自洽计算
3. Electron (Bands+DOS) - 电子能带和态密度

适用场景：
- 新结构的电子性质计算
- 需要先优化再计算能带
- 材料筛选（能带隙、费米面）

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import List, Dict, Any
import logging

from pipelines.base import BasePipeline
from scheduler.task import Task

logger = logging.getLogger(__name__)


class RelaxElectronPipeline(BasePipeline):
    """结构优化 + 电子性质计算Pipeline"""

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
                'name': 'relax',
                'workflow': 'RelaxWorkflow',
                'description': '结构优化',
                'params': {
                    'calculation': 'relax',
                    'nstep': 100,
                },
            },
            {
                'name': 'scf',
                'workflow': 'ScfWorkflow',
                'description': 'SCF自洽计算',
                'params': {
                    'occupations': 'smearing',
                    'smearing': 'mp',
                    'degauss': 0.02,
                },
                'depends_on': ['relax'],
            },
            {
                'name': 'electron_combined',
                'workflow': 'ElectronWorkflow',
                'mode': 'combined',
                'description': '电子能带+DOS计算',
                'params': {
                    'nbnd': None,  # 自动计算
                    'kpoints_dense': (16, 16, 16),  # DOS用密集k点
                },
                'depends_on': ['scf'],
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

            # 收敛参数
            'conv_thr': '1.0d-8',
            'mixing_beta': 0.3,

            # 赝势（示例）
            'pseudopotentials': {},

            # 能带路径（立方系统默认）
            'kpath': [
                [0.0, 0.0, 0.0, 50],  # Gamma
                [0.5, 0.0, 0.5, 50],  # X
                [0.5, 0.5, 0.5, 50],  # M
                [0.0, 0.0, 0.0, 1],   # Gamma
            ],

            # DOS参数
            'dos_emin': -20.0,
            'dos_emax': 20.0,
            'dos_deltae': 0.01,

            # 并行设置
            'nprocs': 8,
        }


def example_usage():
    """使用示例"""
    from pathlib import Path

    # 结构文件列表
    structures = [
        Path("/path/to/structure1.vasp"),
        Path("/path/to/structure2.cif"),
    ]

    # 配置参数
    config = {
        'ecutwfc': 80,
        'ecutrho': 640,
        'kpoints': (8, 8, 8),
        'pseudopotentials': {
            'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
        },
        'nprocs': 8,
    }

    # 创建Pipeline
    pipeline = RelaxElectronPipeline(
        structures=structures,
        work_dir=Path("./calculations"),
        config=config,
    )

    # 运行
    pipeline.run(max_workers=4)


if __name__ == "__main__":
    example_usage()
