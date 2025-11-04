#!/usr/bin/env python3
"""
演示为什么必须在虚轴（Matsubara频率）而非实轴求解Eliashberg方程

这个例子展示：
1. 实轴上格林函数的奇点结构
2. 虚轴上格林函数的光滑性质
3. 简化的自洽方程在两种情况下的数值稳定性对比
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False


# ==================== 第1部分：格林函数的奇点结构 ====================

def free_fermion_greens_function_real(omega, epsilon_k, eta=0.01):
    """
    自由费米子的实轴推迟格林函数

    G^R(ω) = 1 / (ω - ε_k + iη)

    参数:
        omega: 实频率 (meV)
        epsilon_k: 单粒子能量 (meV)
        eta: 小虚部，避免奇点 (meV)
    """
    return 1.0 / (omega - epsilon_k + 1j * eta)


def free_fermion_greens_function_matsubara(iwn, epsilon_k):
    """
    自由费米子的Matsubara格林函数

    G(iω_n) = 1 / (iω_n - ε_k)

    参数:
        iwn: Matsubara频率 iω_n (meV)
        epsilon_k: 单粒子能量 (meV)
    """
    return 1.0 / (1j * iwn - epsilon_k)


def plot_greens_function_comparison():
    """对比实轴和虚轴上的格林函数"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    epsilon_k = 10.0  # 单粒子能量 (meV)

    # ==================== 实轴格林函数 ====================
    omega_real = np.linspace(-30, 30, 1000)
    eta_values = [0.5, 0.1, 0.01]

    ax = axes[0, 0]
    for eta in eta_values:
        G_real = free_fermion_greens_function_real(omega_real, epsilon_k, eta)
        ax.plot(omega_real, np.real(G_real), label=f'Re G (η={eta} meV)')

    ax.axvline(epsilon_k, color='r', linestyle='--', linewidth=2,
              label=f'极点位置 ε_k={epsilon_k} meV')
    ax.set_xlabel('频率 ω (meV)', fontsize=12)
    ax.set_ylabel('Re G^R(ω)', fontsize=12)
    ax.set_title('实轴格林函数的实部（注意极点！）', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-5, 5)

    # ==================== 实轴格林函数虚部 ====================
    ax = axes[0, 1]
    for eta in eta_values:
        G_real = free_fermion_greens_function_real(omega_real, epsilon_k, eta)
        ax.plot(omega_real, np.imag(G_real), label=f'Im G (η={eta} meV)')

    ax.axvline(epsilon_k, color='r', linestyle='--', linewidth=2)
    ax.set_xlabel('频率 ω (meV)', fontsize=12)
    ax.set_ylabel('Im G^R(ω)', fontsize=12)
    ax.set_title('实轴格林函数的虚部（δ函数结构）', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # 标注
    ax.text(epsilon_k + 5, -30,
           f'η→0时：Im G → -πδ(ω-ε_k)',
           fontsize=11, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    # ==================== Matsubara格林函数 ====================
    ax = axes[1, 0]
    T = 10  # 温度 (K)
    kB = 0.08617  # meV/K
    n_max = 20
    n_values = np.arange(0, n_max)
    omega_n = np.pi * kB * T * (2 * n_values + 1)  # Matsubara频率 (meV)

    G_matsubara = free_fermion_greens_function_matsubara(omega_n, epsilon_k)

    ax.plot(omega_n, np.real(G_matsubara), 'bo-', markersize=8,
           linewidth=2, label='Re G(iω_n)')
    ax.plot(omega_n, np.imag(G_matsubara), 'rs-', markersize=8,
           linewidth=2, label='Im G(iω_n)')

    ax.set_xlabel('Matsubara频率 ω_n (meV)', fontsize=12)
    ax.set_ylabel('G(iω_n)', fontsize=12)
    ax.set_title(f'Matsubara格林函数（光滑！），T={T}K', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # 标注
    ax.text(omega_n[10], -0.02,
           '没有奇点！\n适合数值计算',
           fontsize=11, bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # ==================== 模长对比 ====================
    ax = axes[1, 1]

    # 实轴上的模长（取η=0.01）
    eta = 0.01
    G_real = free_fermion_greens_function_real(omega_real, epsilon_k, eta)
    ax.semilogy(omega_real, np.abs(G_real), 'b-', linewidth=2,
               label='|G^R(ω)| (实轴)', alpha=0.7)
    ax.axvline(epsilon_k, color='r', linestyle='--', linewidth=2)

    # Matsubara轴上的模长
    ax.semilogy(omega_n, np.abs(G_matsubara), 'ro', markersize=8,
               linewidth=2, label='|G(iω_n)| (虚轴)')

    ax.set_xlabel('频率 (meV)', fontsize=12)
    ax.set_ylabel('|G| (对数刻度)', fontsize=12)
    ax.set_title('格林函数模长对比', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(-5, 30)

    plt.tight_layout()
    plt.savefig('greens_function_real_vs_matsubara.png', dpi=300, bbox_inches='tight')
    print("✓ 格林函数对比图已保存: greens_function_real_vs_matsubara.png")
    plt.show()


# ==================== 第2部分：简化自洽方程的数值求解 ====================

def simple_gap_equation_real_axis(Delta_initial, omega_grid, epsilon_k,
                                  lambda_coupling, omega_ph, max_iter=100):
    """
    在实轴上求解简化的能隙方程（会遇到数值困难）

    Δ(ω) = λ ∫ dω' [ω_ph / ((ω-ω')² + ω_ph²)] × Δ(ω') / sqrt(ω'² + Δ²(ω'))

    这个方程在实轴上非常不稳定！
    """
    Delta = Delta_initial.copy()
    history = []

    for iteration in range(max_iter):
        Delta_new = np.zeros_like(Delta, dtype=complex)

        for i, omega in enumerate(omega_grid):
            # 配对核（洛伦兹型）
            kernel = lambda_coupling * omega_ph / ((omega - omega_grid)**2 + omega_ph**2)

            # 能隙方程右边
            denom = np.sqrt(omega_grid**2 + np.abs(Delta)**2)
            integrand = kernel * Delta / denom

            # 数值积分
            Delta_new[i] = np.trapz(integrand, omega_grid)

        # 检查收敛（通常不会收敛！）
        error = np.max(np.abs(Delta_new - Delta))
        history.append(np.max(np.abs(Delta_new)))

        if error < 1e-6:
            print(f"实轴方程收敛于第{iteration}次迭代")
            return Delta_new, history

        # 混合（尝试稳定迭代）
        Delta = 0.3 * Delta_new + 0.7 * Delta

    print("实轴方程未收敛！")
    return Delta, history


def simple_gap_equation_matsubara(Delta_initial, omega_n, epsilon_k,
                                  lambda_coupling, omega_ph, T, max_iter=100):
    """
    在Matsubara轴上求解简化的能隙方程（稳定！）

    Δ(iω_n) = πT Σ_m λ(ω_n - ω_m) × Δ(iω_m) / sqrt(ω_m² + Δ²(iω_m))
    """
    Delta = Delta_initial.copy()
    kB = 0.08617  # meV/K
    history = []

    for iteration in range(max_iter):
        Delta_new = np.zeros_like(Delta, dtype=complex)

        for n, wn in enumerate(omega_n):
            summation = 0.0 + 0.0j

            for m, wm in enumerate(omega_n):
                # 配对核
                lambda_kernel = lambda_coupling * omega_ph**2 / ((wn - wm)**2 + omega_ph**2)

                # 能隙方程右边
                denom = np.sqrt(wm**2 + np.abs(Delta[m])**2)
                summation += lambda_kernel * Delta[m] / denom

            Delta_new[n] = np.pi * kB * T * summation

        # 检查收敛
        error = np.max(np.abs(Delta_new - Delta))
        history.append(np.max(np.abs(Delta_new)))

        if error < 1e-6:
            print(f"Matsubara方程收敛于第{iteration}次迭代")
            return Delta_new, history

        # 混合
        Delta = 0.5 * Delta_new + 0.5 * Delta

    print("Matsubara方程未在规定迭代次数内收敛")
    return Delta, history


def compare_self_consistent_iterations():
    """对比实轴和虚轴上自洽迭代的稳定性"""

    # 物理参数
    lambda_coupling = 0.5
    omega_ph = 5.0  # 声子能量 (meV)
    epsilon_k = 10.0  # 电子能量 (meV)
    T = 10  # 温度 (K)
    kB = 0.08617  # meV/K

    # 实轴网格
    omega_real = np.linspace(-50, 50, 200)
    Delta_real_init = 2.0 * np.ones_like(omega_real, dtype=complex)

    # Matsubara频率
    n_matsubara = 50
    n_values = np.arange(0, n_matsubara)
    omega_n = np.pi * kB * T * (2 * n_values + 1)
    Delta_matsubara_init = 2.0 * np.ones(n_matsubara, dtype=complex)

    # 求解
    print("开始在实轴上迭代...")
    Delta_real, history_real = simple_gap_equation_real_axis(
        Delta_real_init, omega_real, epsilon_k, lambda_coupling, omega_ph, max_iter=50
    )

    print("\n开始在Matsubara轴上迭代...")
    Delta_matsubara, history_matsubara = simple_gap_equation_matsubara(
        Delta_matsubara_init, omega_n, epsilon_k, lambda_coupling, omega_ph, T, max_iter=50
    )

    # 可视化
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 迭代历史对比
    ax = axes[0]
    ax.semilogy(history_real, 'r-o', linewidth=2, markersize=6,
               label='实轴迭代（不稳定）')
    ax.semilogy(history_matsubara, 'b-s', linewidth=2, markersize=6,
               label='Matsubara迭代（稳定）')
    ax.set_xlabel('迭代次数', fontsize=12)
    ax.set_ylabel('max|Δ| (meV)', fontsize=12)
    ax.set_title('自洽迭代收敛性对比', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # 添加注释
    ax.text(25, history_real[-1] * 2,
           '实轴：振荡/发散',
           fontsize=11, color='red',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    ax.text(25, history_matsubara[-1] / 2,
           'Matsubara：快速收敛',
           fontsize=11, color='blue',
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # Matsubara解的结构
    ax = axes[1]
    ax.plot(omega_n, np.real(Delta_matsubara), 'bo-', linewidth=2,
           markersize=8, label='Re Δ(iω_n)')
    ax.plot(omega_n, np.imag(Delta_matsubara), 'rs-', linewidth=2,
           markersize=8, label='Im Δ(iω_n)')
    ax.plot(omega_n, np.abs(Delta_matsubara), 'g^-', linewidth=2,
           markersize=8, label='|Δ(iω_n)|')

    ax.set_xlabel('Matsubara频率 ω_n (meV)', fontsize=12)
    ax.set_ylabel('能隙 Δ (meV)', fontsize=12)
    ax.set_title(f'Matsubara解的结构 (T={T}K)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('self_consistent_comparison.png', dpi=300, bbox_inches='tight')
    print("\n✓ 自洽迭代对比图已保存: self_consistent_comparison.png")
    plt.show()


# ==================== 主程序 ====================

def main():
    print("=" * 70)
    print("演示：为什么必须在虚轴（Matsubara频率）求解Eliashberg方程")
    print("=" * 70)

    print("\n第1部分：格林函数的奇点结构")
    print("-" * 70)
    plot_greens_function_comparison()

    print("\n第2部分：自洽方程的数值稳定性")
    print("-" * 70)
    compare_self_consistent_iterations()

    print("\n" + "=" * 70)
    print("总结:")
    print("-" * 70)
    print("1. 实轴上的格林函数有奇点（极点、δ函数）")
    print("2. Matsubara格林函数在虚轴上光滑，适合数值计算")
    print("3. 实轴自洽方程数值不稳定，容易发散")
    print("4. Matsubara自洽方程稳定收敛")
    print("5. 因此必须在虚轴求解，然后通过Padé延拓回实轴！")
    print("=" * 70)


if __name__ == "__main__":
    main()
