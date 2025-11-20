"""
超导Tc计算模块

实现McMillan、Allen-Dynes、Eliashberg方程的Tc计算。

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import Tuple, Optional
import numpy as np
from core.constants import BOLTZMANN_CONSTANT, ELECTRON_CHARGE


class TcCalculator:
    """超导临界温度计算器"""

    @staticmethod
    def calculate_lambda(frequencies: np.ndarray, alpha2f: np.ndarray) -> Tuple[float, float]:
        """
        计算电声耦合常数λ和对数平均声子频率ω_log

        Parameters
        ----------
        frequencies : np.ndarray
            频率数组 [Ry]
        alpha2f : np.ndarray
            Alpha2F函数

        Returns
        -------
        lambda_value : float
            电声耦合常数
        omega_log : float
            对数平均声子频率 [K]
        """
        # λ = 2 ∫ α²F(ω)/ω dω
        integrand_lambda = alpha2f / (frequencies + 1e-10)  # 避免除零
        lambda_value = 2.0 * np.trapz(integrand_lambda, frequencies)

        # ω_log = exp[2/λ ∫ (α²F(ω)/ω) ln(ω) dω]
        integrand_omega = (alpha2f / (frequencies + 1e-10)) * np.log(frequencies + 1e-10)
        omega_log_ry = np.exp((2.0 / lambda_value) * np.trapz(integrand_omega, frequencies))

        # 转换为K（1 Ry = 157887 K）
        omega_log = omega_log_ry * 157887.0

        return lambda_value, omega_log

    @staticmethod
    def mcmillan_tc(lambda_value: float, omega_log: float, mu: float = 0.10) -> float:
        """
        McMillan公式计算Tc

        Tc = (ω_log / 1.2) * exp[-1.04(1+λ) / (λ - μ*(1+0.62λ))]

        Parameters
        ----------
        lambda_value : float
            电声耦合常数
        omega_log : float
            对数平均声子频率 [K]
        mu : float
            屏蔽常数 μ*

        Returns
        -------
        float
            Tc [K]
        """
        numerator = 1.04 * (1.0 + lambda_value)
        denominator = lambda_value - mu * (1.0 + 0.62 * lambda_value)

        if denominator <= 0:
            return 0.0

        tc = (omega_log / 1.2) * np.exp(-numerator / denominator)
        return tc

    @staticmethod
    def allen_dynes_tc(lambda_value: float, omega_log: float, mu: float = 0.10) -> float:
        """
        Allen-Dynes公式计算Tc

        考虑了强耦合和Coulomb赝势修正

        Parameters
        ----------
        lambda_value : float
            电声耦合常数
        omega_log : float
            对数平均声子频率 [K]
        mu : float
            屏蔽常数 μ*

        Returns
        -------
        float
            Tc [K]
        """
        # f1 和 f2 修正因子
        f1 = (1.0 + (lambda_value / 2.46) ** (3.0/2.0)) ** (1.0/3.0)
        f2 = 1.0 + (omega_log / 250.0 - 1.0) * (lambda_value / (lambda_value + 1.0))

        # 有效耦合常数
        lambda_eff = lambda_value - mu * (1.0 + 0.62 * lambda_value)

        if lambda_eff <= 0:
            return 0.0

        # Tc计算
        tc = (f1 * f2 * omega_log / 1.2) * np.exp(-1.04 * (1.0 + lambda_value) / lambda_eff)

        return tc

    @staticmethod
    def read_alpha2f(alpha2f_file: Path) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取ALPHA2F.OUT文件

        Parameters
        ----------
        alpha2f_file : Path
            ALPHA2F.OUT文件路径

        Returns
        -------
        frequencies : np.ndarray
            频率 [Hartree]
        alpha2f : np.ndarray
            Alpha2F函数
        """
        data = np.loadtxt(alpha2f_file)
        frequencies = data[:, 0]  # Hartree
        alpha2f = data[:, 1]

        # 转换为Rydberg（QE内部单位）
        frequencies_ry = frequencies * 2.0  # 1 Ha = 2 Ry

        return frequencies_ry, alpha2f


# 使用示例
if __name__ == "__main__":
    calc = TcCalculator()

    # 示例：计算H3S的Tc
    # frequencies = ...  # 从文件读取
    # alpha2f = ...
    # lambda_val, omega_log = calc.calculate_lambda(frequencies, alpha2f)
    # tc_mcm = calc.mcmillan_tc(lambda_val, omega_log, mu=0.10)
    # tc_ad = calc.allen_dynes_tc(lambda_val, omega_log, mu=0.10)
    # print(f"λ = {lambda_val:.3f}")
    # print(f"ω_log = {omega_log:.1f} K")
    # print(f"Tc (McMillan) = {tc_mcm:.1f} K")
    # print(f"Tc (Allen-Dynes) = {tc_ad:.1f} K")
    pass
