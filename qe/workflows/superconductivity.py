"""
超导性质计算工作流

支持超导计算的完整流程：
- lambda: 计算电声耦合常数λ
- tc: 计算超导临界温度Tc（McMillan/Allen-Dynes/Eliashberg）
- alpha2f: 提取并分析Eliashberg谱函数α²F(ω)

作者：Claude
创建时间：2025-11-20
"""

from pathlib import Path
from typing import Optional, List, Dict, Tuple
import logging
import numpy as np

from workflows.base import BaseWorkflow
from config.superconductivity import SuperconductivityConfig
from analysis.superconductivity.tc_calculator import TcCalculator

logger = logging.getLogger(__name__)


class SuperconductivityWorkflow(BaseWorkflow):
    """超导性质计算工作流"""

    def __init__(self, config: SuperconductivityConfig, mode: str = 'lambda'):
        """
        初始化超导工作流

        Parameters
        ----------
        config : SuperconductivityConfig
            超导配置对象
        mode : str
            计算模式: lambda, tc, alpha2f
        """
        super().__init__(config)
        self.mode = mode
        self.tc_calculator = TcCalculator()

    def prepare_input(self) -> Path:
        """根据mode准备相应的输入文件"""
        if self.mode == 'lambda':
            return self._prepare_lambda_input()
        elif self.mode == 'tc':
            return self._prepare_tc_input()
        elif self.mode == 'alpha2f':
            return self._prepare_alpha2f_input()
        else:
            raise ValueError(f"未知的超导计算模式: {self.mode}")

    def _prepare_lambda_input(self) -> Path:
        """准备lambda.x计算输入"""
        input_file = self.work_dir / 'lambda.in'

        with open(input_file, 'w') as f:
            f.write("&INPUT\n")
            f.write(f"  prefix = '{self.work_dir.name}'\n")
            f.write("  outdir = './tmp/'\n")
            f.write("  fildyn = 'dyn'\n")

            # 电声耦合参数
            sigma = self.config.get_param('el_ph_sigma', 0.005)
            nsigma = self.config.get_param('el_ph_nsigma', 10)

            f.write(f"  el_ph_sigma = {sigma}\n")
            f.write(f"  el_ph_nsigma = {nsigma}\n")

            # 是否计算Gamma线宽
            if self.config.get_param('compute_gamma'):
                f.write("  la2F = .true.\n")

            f.write("/\n")

        logger.info(f"Lambda计算输入已生成: {input_file}")
        return input_file

    def _prepare_tc_input(self) -> Path:
        """
        准备Tc计算输入（实际是Python脚本，使用已有的alpha2f数据）

        Returns
        -------
        Path
            Tc计算脚本路径
        """
        script_file = self.work_dir / 'compute_tc.py'

        with open(script_file, 'w') as f:
            f.write('#!/usr/bin/env python3\n')
            f.write('"""\n')
            f.write('超导Tc计算脚本\n')
            f.write('读取alpha2F数据并使用McMillan和Allen-Dynes公式计算Tc\n')
            f.write('"""\n\n')

            f.write('import numpy as np\n')
            f.write('from pathlib import Path\n')
            f.write('import sys\n\n')

            f.write('# 添加项目路径\n')
            f.write('sys.path.append(str(Path(__file__).parent.parent.parent))\n\n')

            f.write('from analysis.superconductivity.tc_calculator import TcCalculator\n\n')

            f.write('def main():\n')
            f.write('    calc = TcCalculator()\n\n')

            f.write('    # 读取alpha2F文件\n')
            f.write('    alpha2f_file = Path("_ph0") / "fildyn.a2Fsave"\n')
            f.write('    if not alpha2f_file.exists():\n')
            f.write('        print(f"错误: 找不到alpha2F文件 {alpha2f_file}")\n')
            f.write('        return\n\n')

            f.write('    frequencies, alpha2f = calc.read_alpha2f(alpha2f_file)\n\n')

            f.write('    # 计算λ和ω_log\n')
            f.write('    lambda_val, omega_log = calc.calculate_lambda(frequencies, alpha2f)\n\n')

            f.write('    print(f"电声耦合常数 λ = {lambda_val:.4f}")\n')
            f.write('    print(f"对数平均声子频率 ω_log = {omega_log:.2f} K")\n')
            f.write('    print()\n\n')

            # 获取mu值列表
            mu_values = self.config.get_param('mu_values', [0.10, 0.13])
            f.write(f'    mu_values = {mu_values}\n\n')

            f.write('    # 使用不同公式和mu值计算Tc\n')
            f.write('    results = []\n')
            f.write('    for mu in mu_values:\n')
            f.write('        tc_mcm = calc.mcmillan_tc(lambda_val, omega_log, mu)\n')
            f.write('        tc_ad = calc.allen_dynes_tc(lambda_val, omega_log, mu)\n\n')

            f.write('        results.append({\n')
            f.write('            "mu": mu,\n')
            f.write('            "tc_mcmillan": tc_mcm,\n')
            f.write('            "tc_allen_dynes": tc_ad\n')
            f.write('        })\n\n')

            f.write('        print(f"μ* = {mu:.2f}:")\n')
            f.write('        print(f"  Tc (McMillan)    = {tc_mcm:6.2f} K")\n')
            f.write('        print(f"  Tc (Allen-Dynes) = {tc_ad:6.2f} K")\n')
            f.write('        print()\n\n')

            f.write('    # 保存结果到文件\n')
            f.write('    with open("tc_results.txt", "w") as out:\n')
            f.write('        out.write(f"电声耦合常数 λ = {lambda_val:.4f}\\n")\n')
            f.write('        out.write(f"对数平均声子频率 ω_log = {omega_log:.2f} K\\n\\n")\n')
            f.write('        out.write("超导临界温度 Tc:\\n")\n')
            f.write('        out.write("-" * 60 + "\\n")\n')
            f.write('        out.write(f"{"μ*":>6s} {"McMillan (K)":>15s} {"Allen-Dynes (K)":>18s}\\n")\n')
            f.write('        out.write("-" * 60 + "\\n")\n\n')

            f.write('        for res in results:\n')
            f.write('            out.write(f\'{res["mu"]:6.2f} {res["tc_mcmillan"]:15.2f} {res["tc_allen_dynes"]:18.2f}\\n\')\n\n')

            f.write('    print("结果已保存到 tc_results.txt")\n\n')

            f.write('if __name__ == "__main__":\n')
            f.write('    main()\n')

        script_file.chmod(0o755)
        logger.info(f"Tc计算脚本已生成: {script_file}")
        return script_file

    def _prepare_alpha2f_input(self) -> Path:
        """准备提取alpha2F数据的脚本"""
        script_file = self.work_dir / 'extract_alpha2f.py'

        with open(script_file, 'w') as f:
            f.write('#!/usr/bin/env python3\n')
            f.write('"""\n')
            f.write('提取和分析Eliashberg谱函数α²F(ω)\n')
            f.write('"""\n\n')

            f.write('import numpy as np\n')
            f.write('import matplotlib.pyplot as plt\n')
            f.write('from pathlib import Path\n\n')

            f.write('def read_alpha2f(filename):\n')
            f.write('    """读取alpha2F文件"""\n')
            f.write('    data = np.loadtxt(filename)\n')
            f.write('    omega = data[:, 0]  # 频率 (Hartree)\n')
            f.write('    alpha2f = data[:, 1]  # α²F(ω)\n')
            f.write('    return omega, alpha2f\n\n')

            f.write('def plot_alpha2f(omega, alpha2f, output="alpha2f.png"):\n')
            f.write('    """绘制α²F(ω)图"""\n')
            f.write('    # 转换单位：Hartree -> meV\n')
            f.write('    omega_mev = omega * 27211.386  # 1 Ha = 27211.386 meV\n\n')

            f.write('    plt.figure(figsize=(8, 6), dpi=300)\n')
            f.write('    plt.plot(omega_mev, alpha2f, linewidth=2)\n')
            f.write('    plt.xlabel("Frequency ω (meV)", fontsize=14)\n')
            f.write('    plt.ylabel("α²F(ω)", fontsize=14)\n')
            f.write('    plt.title("Eliashberg Spectral Function", fontsize=16)\n')
            f.write('    plt.grid(alpha=0.3)\n')
            f.write('    plt.tight_layout()\n')
            f.write('    plt.savefig(output)\n')
            f.write('    print(f"α²F(ω)图已保存: {output}")\n\n')

            f.write('def main():\n')
            f.write('    alpha2f_file = Path("_ph0") / "fildyn.a2Fsave"\n\n')

            f.write('    if not alpha2f_file.exists():\n')
            f.write('        print(f"错误: 找不到alpha2F文件 {alpha2f_file}")\n')
            f.write('        return\n\n')

            f.write('    omega, alpha2f = read_alpha2f(alpha2f_file)\n\n')

            f.write('    # 保存数据\n')
            f.write('    np.savetxt("alpha2f_data.txt", np.column_stack([omega, alpha2f]),\n')
            f.write('               header="Frequency(Ha)  alpha2F", fmt="%15.8e")\n\n')

            f.write('    # 绘图\n')
            f.write('    plot_alpha2f(omega, alpha2f)\n\n')

            f.write('    print(f"频率范围: {omega.min():.6f} - {omega.max():.6f} Ha")\n')
            f.write('    print(f"α²F最大值: {alpha2f.max():.6f}")\n\n')

            f.write('if __name__ == "__main__":\n')
            f.write('    main()\n')

        script_file.chmod(0o755)
        logger.info(f"Alpha2F提取脚本已生成: {script_file}")
        return script_file

    def generate_submit_script(self) -> Path:
        """生成超导计算提交脚本"""
        if self.mode == 'lambda':
            return self._generate_lambda_script()
        elif self.mode == 'tc':
            return self._generate_tc_script()
        elif self.mode == 'alpha2f':
            return self._generate_alpha2f_script()
        else:
            raise ValueError(f"未知的超导计算模式: {self.mode}")

    def _generate_lambda_script(self) -> Path:
        """生成lambda计算脚本"""
        script_file = self.work_dir / 's_lambda.sh'

        nprocs = self.config.get_param('nprocs', 8)

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 电声耦合常数λ计算\n\n")

            f.write("echo '计算电声耦合常数λ'\n")
            f.write(f"mpirun -np {nprocs} lambda.x < lambda.in > lambda.out\n\n")

            f.write("# 检查输出文件\n")
            f.write("if [ -f lambda.out ]; then\n")
            f.write('    lambda_val=$(grep "lambda" lambda.out | tail -1 | awk \'{print $NF}\')\n')
            f.write('    echo "电声耦合常数 λ = $lambda_val"\n')
            f.write("else\n")
            f.write("    echo '错误: lambda计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            # 检查是否生成了alpha2F文件
            f.write("if [ -f _ph0/*.a2Fsave ]; then\n")
            f.write("    echo 'α²F(ω)数据已生成'\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"Lambda计算脚本已生成: {script_file}")
        return script_file

    def _generate_tc_script(self) -> Path:
        """生成Tc计算脚本"""
        script_file = self.work_dir / 's_tc.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# 超导临界温度Tc计算\n\n")

            f.write("echo '计算超导临界温度Tc'\n\n")

            # 检查alpha2F文件是否存在
            f.write("if [ ! -f _ph0/fildyn.a2Fsave ]; then\n")
            f.write("    echo '错误: 找不到alpha2F文件！请先运行lambda计算'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("# 运行Tc计算脚本\n")
            f.write("python3 compute_tc.py\n\n")

            f.write("if [ -f tc_results.txt ]; then\n")
            f.write("    echo ''\n")
            f.write("    echo '========== Tc计算结果 =========='\n")
            f.write("    cat tc_results.txt\n")
            f.write("else\n")
            f.write("    echo '错误: Tc计算失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"Tc计算脚本已生成: {script_file}")
        return script_file

    def _generate_alpha2f_script(self) -> Path:
        """生成alpha2F提取和绘图脚本"""
        script_file = self.work_dir / 's_alpha2f.sh'

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# α²F(ω)提取和分析\n\n")

            f.write("echo '提取α²F(ω)数据'\n\n")

            # 检查alpha2F文件
            f.write("if [ ! -f _ph0/fildyn.a2Fsave ]; then\n")
            f.write("    echo '错误: 找不到alpha2F文件！'\n")
            f.write("    exit 1\n")
            f.write("fi\n\n")

            f.write("# 运行提取脚本\n")
            f.write("python3 extract_alpha2f.py\n\n")

            f.write("if [ -f alpha2f_data.txt ] && [ -f alpha2f.png ]; then\n")
            f.write("    echo 'α²F(ω)数据和图像已生成'\n")
            f.write("else\n")
            f.write("    echo '错误: α²F提取失败！'\n")
            f.write("    exit 1\n")
            f.write("fi\n")

        script_file.chmod(0o755)
        logger.info(f"Alpha2F提取脚本已生成: {script_file}")
        return script_file

    def calculate_tc_from_data(self, alpha2f_file: Path, mu_values: Optional[List[float]] = None) -> Dict:
        """
        从alpha2F文件直接计算Tc

        Parameters
        ----------
        alpha2f_file : Path
            alpha2F数据文件路径
        mu_values : List[float], optional
            屏蔽常数μ*列表

        Returns
        -------
        Dict
            包含lambda, omega_log和各mu值对应的Tc结果
        """
        if mu_values is None:
            mu_values = self.config.get_param('mu_values', [0.10, 0.13])

        # 读取alpha2F数据
        frequencies, alpha2f = self.tc_calculator.read_alpha2f(alpha2f_file)

        # 计算λ和ω_log
        lambda_val, omega_log = self.tc_calculator.calculate_lambda(frequencies, alpha2f)

        results = {
            'lambda': lambda_val,
            'omega_log': omega_log,
            'tc_values': []
        }

        # 计算不同mu值的Tc
        for mu in mu_values:
            tc_mcm = self.tc_calculator.mcmillan_tc(lambda_val, omega_log, mu)
            tc_ad = self.tc_calculator.allen_dynes_tc(lambda_val, omega_log, mu)

            results['tc_values'].append({
                'mu': mu,
                'tc_mcmillan': tc_mcm,
                'tc_allen_dynes': tc_ad
            })

        return results
