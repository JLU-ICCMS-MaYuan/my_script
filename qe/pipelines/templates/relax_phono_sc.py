"""
结构优化→声子→超导 Pipeline模板

作者：Claude
创建时间：2025-11-19
"""

from pipelines.base import BasePipeline


class RelaxPhonoSuperconductivityPipeline(BasePipeline):
    """
    完整的超导性质计算流程

    流程步骤：
    1. 结构优化 (relax-vc)
    2. 密K点自洽 (scffit)
    3. 稀K点自洽 (scf)
    4. 声子计算 (phonon)
    5. 声子态密度 (phonodos)
    6. Lambda计算
    7. Tc计算 (McMillan/Allen-Dynes/Eliashberg)
    """

    def define_steps(self):
        """定义计算步骤"""
        qpoints = self.config.get('qpoints', '6 6 6')
        ecutwfc = self.config.get('ecutwfc', 80)
        mu_values = self.config.get('mu_values', [0.10, 0.13])

        steps = [
            {
                'name': 'relax',
                'workflow': 'RelaxWorkflow',
                'params': {
                    'mode': 'relax-vc',
                    'ecutwfc': ecutwfc,
                    'kpoints': qpoints
                }
            },
            {
                'name': 'scffit',
                'workflow': 'ScfWorkflow',
                'params': {
                    'mode': 'scffit',
                    'ecutwfc': ecutwfc,
                    'kpoints_dense': f"{int(qpoints.split()[0])*4} {int(qpoints.split()[1])*4} {int(qpoints.split()[2])*4}"
                }
            },
            {
                'name': 'scf',
                'workflow': 'ScfWorkflow',
                'params': {
                    'mode': 'scf',
                    'ecutwfc': ecutwfc,
                    'kpoints': qpoints
                }
            },
            {
                'name': 'phonon',
                'workflow': 'PhononWorkflow',
                'params': {
                    'mode': 'nosplit',
                    'qpoints': qpoints
                }
            },
            {
                'name': 'phonodos',
                'workflow': 'PhononWorkflow',
                'params': {
                    'mode': 'phonodos',
                    'qpoints_dense': f"{int(qpoints.split()[0])*2} {int(qpoints.split()[1])*2} {int(qpoints.split()[2])*2}"
                }
            },
            {
                'name': 'lambda',
                'workflow': 'SuperconductivityWorkflow',
                'params': {
                    'mode': 'lambda'
                }
            },
            {
                'name': 'tc',
                'workflow': 'SuperconductivityWorkflow',
                'params': {
                    'mode': 'tc',
                    'mu_values': mu_values
                }
            }
        ]

        return steps
