#!/usr/bin/env python3
"""
测试Allen-Dynes两种f2形式在理想情况下的一致性
当ω_rms/ω_log接近1时，两种形式应该给出几乎相同的结果
"""

import numpy as np

def test_allen_dynes_convergence():
    """测试当频率比接近1时两种形式的收敛性"""

    # 测试参数
    lam = 1.0
    mustar = 0.16

    print("Allen-Dynes f2因子收敛性测试")
    print("="*50)
    print(f"λ = {lam:.3f}, μ* = {mustar:.3f}")
    print()
    print("ω_rms/ω_log    Fig1_f2    Fig2_f2    Tc1(K)    Tc2(K)    相对误差(%)")
    print("-" * 70)

    # 基础临界温度参数 (假设ω_log = 0.003 Ha)
    omega_log = 0.003
    tc_base = omega_log / (1.2 * 3.166815343e-6) * np.exp(-1.04 * (1 + lam) / (lam - mustar - 0.62 * lam * mustar))

    # f1因子 (两种形式相同)
    l1 = 2.46 * (1.0 + 3.8 * mustar)
    f1 = (1.0 + (lam / l1) ** 1.5) ** (1.0 / 3.0)

    # 测试不同的频率比
    ratios = [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5]

    for omega_ratio in ratios:
        # Fig1形式
        l2_fig1 = 1.82 * (1.0 + 6.3 * mustar) * omega_ratio
        f2_fig1 = 1.0 + (omega_ratio - 1.0) * lam**2 / (lam**2 + l2_fig1**2)
        tc1 = tc_base * f1 * f2_fig1

        # Fig2形式
        l2_fig2_squared = 3.312 * (1.0 + 6.3 * mustar)**2
        f2_fig2 = 1.0 - lam**2 * (1.0 - omega_ratio) / (lam**2 + l2_fig2_squared)
        tc2 = tc_base * f1 * f2_fig2

        # 相对误差
        rel_error = abs(tc1 - tc2) / tc1 * 100

        print(f"{omega_ratio:8.2f}    {f2_fig1:8.5f}    {f2_fig2:8.5f}    {tc1:8.2f}    {tc2:8.2f}    {rel_error:8.4f}")

    print()
    print("结论:")
    print("- 当ω_rms/ω_log = 1.00时，两种形式完全一致")
    print("- 偏离1越多，误差越大")
    print("- 误差<1%时可认为两种形式等价")
    print("- 误差>2%时建议使用Fig1精确形式")

if __name__ == "__main__":
    test_allen_dynes_convergence()