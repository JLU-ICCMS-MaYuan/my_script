#!/usr/bin/env python3
"""
Padé解析延拓的完整实现与验证

演示内容：
1. Vidberg-Serene递归算法的数学原理
2. 从Matsubara频率到实轴的延拓过程
3. 用已知解析函数验证Padé近似的准确性
4. 超导能隙函数的实际延拓示例
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False


# ==================== 第1部分：Padé近似的基础实现 ====================

def vidberg_serene_pade(z_in, f_in, z_out):
    """
    Vidberg-Serene递归算法实现Padé解析延拓

    给定函数在一组点{z_in[i], f_in[i]}上的值，
    通过Padé近似计算函数在z_out处的值

    参数:
        z_in: 输入点（复数数组），如Matsubara频率 iω_n
        f_in: 函数在输入点的值
        z_out: 输出点（复数或实数数组）

    返回:
        f_out: Padé近似在z_out处的值

    数学原理：
        构造有理函数 P(z) = p(z)/q(z) 使得 P(z_in[i]) = f_in[i]
        使用递归公式：
            g[0,j] = f_in[j]
            g[i,j] = (g[i-1,i-1] - g[i-1,j]) / [(z_in[j] - z_in[i-1]) * g[i-1,j]]
    """
    n_in = len(z_in)
    n_out = len(z_out) if isinstance(z_out, np.ndarray) else 1

    # 如果z_out是标量，转换为数组
    if not isinstance(z_out, np.ndarray):
        z_out = np.array([z_out])

    # 初始化输出
    f_out = np.zeros(n_out, dtype=complex)

    # 构建g矩阵（递归计算）
    g = np.zeros((n_in, n_in), dtype=complex)
    g[0, :] = f_in  # 初始化第一行

    # Vidberg-Serene递归
    for i in range(1, n_in):
        for j in range(i, n_in):
            numerator = g[i-1, i-1] - g[i-1, j]
            denominator = (z_in[j] - z_in[i-1]) * g[i-1, j]

            # 避免除零
            if np.abs(denominator) < 1e-15:
                g[i, j] = 0.0
            else:
                g[i, j] = numerator / denominator

    # 对每个输出点计算Padé近似（连分式展开）
    for idx, z in enumerate(z_out):
        # 初始化连分式递归
        a0, a1 = 0.0 + 0.0j, g[0, 0]
        b0, b1 = 1.0 + 0.0j, 1.0 + 0.0j

        # 递归计算连分式
        for j in range(1, n_in):
            z_diff_times_g = (z - z_in[j-1]) * g[j, j]

            # 连分式递归关系
            a_new = a1 + z_diff_times_g * a0
            b_new = b1 + z_diff_times_g * b0

            # 更新
            a0, a1 = a1, a_new
            b0, b1 = b1, b_new

        # Padé近似 = a1 / b1
        if np.abs(b1) > 1e-15:
            f_out[idx] = a1 / b1
        else:
            f_out[idx] = np.nan

    return f_out[0] if n_out == 1 else f_out


# ==================== 第2部分：用已知函数验证Padé近似 ====================

def test_pade_with_known_function():
    """
    用已知解析函数测试Padé延拓的准确性

    测试函数：BCS型能隙函数
        Δ(z) = Δ₀ * z / sqrt(z² - Δ₀²)

    在Matsubara点采样，然后延拓到实轴对比
    """
    print("\n" + "="*70)
    print("第1部分：用已知函数验证Padé近似")
    print("="*70)

    # 参数
    Delta_0 = 10.0  # 超导能隙 (meV)
    T = 10.0  # 温度 (K)
    kB = 0.08617  # meV/K

    # BCS能隙函数（已知解析形式）
    def delta_bcs(z):
        """BCS能隙函数"""
        return Delta_0 * z / np.sqrt(z**2 - Delta_0**2 + 0j)

    # 生成Matsubara频率点
    n_matsubara = 30
    n_values = np.arange(0, n_matsubara)
    omega_n = np.pi * kB * T * (2 * n_values + 1)  # 正频率
    z_matsubara = 1j * omega_n  # iω_n

    # 在Matsubara点计算精确值
    delta_matsubara = delta_bcs(z_matsubara)

    print(f"✓ 在{n_matsubara}个Matsubara点采样")
    print(f"  温度 T = {T} K")
    print(f"  能隙 Δ₀ = {Delta_0} meV")

    # 延拓到实轴
    omega_real = np.linspace(0.1, 50, 200)
    z_real = omega_real + 1j * 0.01  # ω + i0⁺

    # Padé延拓
    delta_pade = vidberg_serene_pade(z_matsubara, delta_matsubara, z_real)

    # 精确解
    delta_exact = delta_bcs(z_real)

    # 计算误差
    relative_error = np.abs((delta_pade - delta_exact) / delta_exact)

    print(f"✓ Padé延拓完成")
    print(f"  平均相对误差: {np.mean(relative_error):.6f}")
    print(f"  最大相对误差: {np.max(relative_error):.6f}")

    # 可视化
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 图1：Matsubara轴上的采样点
    ax = axes[0, 0]
    ax.plot(np.imag(z_matsubara), np.real(delta_matsubara), 'bo-',
           markersize=8, linewidth=2, label='Re Δ(iω_n)')
    ax.plot(np.imag(z_matsubara), np.imag(delta_matsubara), 'rs-',
           markersize=8, linewidth=2, label='Im Δ(iω_n)')
    ax.set_xlabel('Matsubara频率 ω_n (meV)', fontsize=12)
    ax.set_ylabel('Δ(iω_n) (meV)', fontsize=12)
    ax.set_title('Matsubara轴上的采样', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # 图2：实轴延拓结果对比
    ax = axes[0, 1]
    ax.plot(omega_real, np.real(delta_exact), 'b-',
           linewidth=3, label='精确解 Re Δ(ω)', alpha=0.5)
    ax.plot(omega_real, np.real(delta_pade), 'r--',
           linewidth=2, label='Padé近似 Re Δ(ω)')
    ax.axhline(Delta_0, color='green', linestyle='--', linewidth=1.5,
              label=f'Δ₀ = {Delta_0} meV')
    ax.set_xlabel('实频率 ω (meV)', fontsize=12)
    ax.set_ylabel('Re Δ(ω) (meV)', fontsize=12)
    ax.set_title('实轴延拓：实部对比', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-5, 15)

    # 图3：虚部对比
    ax = axes[1, 0]
    ax.plot(omega_real, np.imag(delta_exact), 'b-',
           linewidth=3, label='精确解 Im Δ(ω)', alpha=0.5)
    ax.plot(omega_real, np.imag(delta_pade), 'r--',
           linewidth=2, label='Padé近似 Im Δ(ω)')
    ax.set_xlabel('实频率 ω (meV)', fontsize=12)
    ax.set_ylabel('Im Δ(ω) (meV)', fontsize=12)
    ax.set_title('实轴延拓：虚部对比', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # 图4：相对误差
    ax = axes[1, 1]
    ax.semilogy(omega_real, relative_error, 'g-', linewidth=2)
    ax.set_xlabel('实频率 ω (meV)', fontsize=12)
    ax.set_ylabel('相对误差', fontsize=12)
    ax.set_title('Padé近似的精度', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')

    # 标注
    ax.text(25, np.max(relative_error) / 2,
           f'平均误差: {np.mean(relative_error):.2e}\n最大误差: {np.max(relative_error):.2e}',
           fontsize=11, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    plt.tight_layout()
    plt.savefig('pade_validation_known_function.png', dpi=300, bbox_inches='tight')
    print("✓ 图像已保存: pade_validation_known_function.png")
    plt.show()


# ==================== 第3部分：从Eliashberg解延拓到实轴 ====================

def example_eliashberg_continuation():
    """
    演示从Eliashberg方程的Matsubara解延拓到实轴

    这是eliashberg_solver.py中实际使用的场景
    """
    print("\n" + "="*70)
    print("第2部分：Eliashberg能隙的实际延拓")
    print("="*70)

    # 模拟Eliashberg求解结果（Matsubara轴上的解）
    T = 50.0  # K
    kB = 0.08617  # meV/K
    Delta_0 = 15.0  # meV

    # Matsubara频率
    n_max = 40
    n_values = np.arange(0, n_max)
    omega_n = np.pi * kB * T * (2 * n_values + 1)

    # 模拟典型的Eliashberg解（BCS形式 + 频率依赖修正）
    Z_n = 1.0 + 0.3 * Delta_0 / omega_n  # 重整化因子
    Delta_n = Delta_0 * omega_n / np.sqrt(omega_n**2 + (Delta_0/Z_n)**2)

    print(f"✓ 模拟Matsubara解")
    print(f"  温度 T = {T} K")
    print(f"  能隙幅度 Δ₀ ≈ {Delta_0} meV")
    print(f"  Matsubara点数: {n_max}")

    # 延拓到实轴
    omega_real = np.linspace(0.1, 100, 300)
    z_real = omega_real + 1j * 0.1

    # 分别延拓Δ和Z
    z_matsubara = 1j * omega_n
    Delta_real = vidberg_serene_pade(z_matsubara, Delta_n, z_real)
    Z_real = vidberg_serene_pade(z_matsubara, Z_n, z_real)

    # 计算谱函数 A(ω) = Im[Z(ω)] / [ω² + Re[Δ(ω)]²]^(1/2)
    # 这是ARPES等实验观测的量
    spectral_function = np.abs(np.imag(Z_real)) / np.sqrt(
        omega_real**2 + np.real(Delta_real)**2
    )

    print(f"✓ 延拓完成")
    print(f"  实频点数: {len(omega_real)}")

    # 可视化
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 图1：能隙函数Δ
    ax = axes[0, 0]
    ax.plot(omega_n, np.real(Delta_n), 'bo', markersize=6,
           label='Matsubara: Re Δ(iω_n)')
    ax.plot(omega_real, np.real(Delta_real), 'r-', linewidth=2,
           label='实轴延拓: Re Δ(ω)', alpha=0.7)
    ax.axhline(Delta_0, color='green', linestyle='--',
              label=f'Δ₀ ≈ {Delta_0} meV')
    ax.set_xlabel('频率 (meV)', fontsize=12)
    ax.set_ylabel('Re Δ (meV)', fontsize=12)
    ax.set_title('超导能隙函数', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 100)

    # 图2：重整化因子Z
    ax = axes[0, 1]
    ax.plot(omega_n, np.real(Z_n), 'bs', markersize=6,
           label='Matsubara: Re Z(iω_n)')
    ax.plot(omega_real, np.real(Z_real), 'r-', linewidth=2,
           label='实轴延拓: Re Z(ω)', alpha=0.7)
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('频率 (meV)', fontsize=12)
    ax.set_ylabel('Re Z', fontsize=12)
    ax.set_title('准粒子重整化因子', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 100)

    # 图3：谱函数A(ω)
    ax = axes[1, 0]
    ax.plot(omega_real, spectral_function, 'b-', linewidth=2)
    ax.fill_between(omega_real, spectral_function, alpha=0.3)
    ax.set_xlabel('频率 ω (meV)', fontsize=12)
    ax.set_ylabel('谱函数 A(ω) (任意单位)', fontsize=12)
    ax.set_title('准粒子谱函数（ARPES可观测）', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 50)

    # 标注能隙边缘
    gap_edge = Delta_0
    ax.axvline(gap_edge, color='r', linestyle='--', linewidth=2,
              label=f'能隙边缘 ≈ {gap_edge:.1f} meV')
    ax.legend(fontsize=10)

    # 图4：复平面上的延拓路径
    ax = axes[1, 1]

    # Matsubara点
    ax.plot([0]*len(omega_n), omega_n, 'bo', markersize=8,
           label='Matsubara点 iω_n', zorder=3)

    # 实轴目标点
    ax.plot(omega_real, [0.1]*len(omega_real), 'r.', markersize=2,
           label='实轴目标 ω+i0⁺', alpha=0.5)

    # 延拓路径示意
    for i in range(0, n_max, 5):
        path_re = [0, omega_n[i]*0.8]
        path_im = [omega_n[i], 0.1]
        ax.plot(path_re, path_im, 'g--', alpha=0.2, linewidth=1)

    # 坐标轴
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    ax.set_xlabel('Re z (meV)', fontsize=12)
    ax.set_ylabel('Im z (meV)', fontsize=12)
    ax.set_title('复平面上的解析延拓', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-10, 100)
    ax.set_ylim(-5, 50)

    # 标注
    ax.text(5, 40, 'Padé延拓\n从虚轴→实轴', fontsize=11,
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    plt.tight_layout()
    plt.savefig('pade_eliashberg_continuation.png', dpi=300, bbox_inches='tight')
    print("✓ 图像已保存: pade_eliashberg_continuation.png")
    plt.show()


# ==================== 第4部分：Padé近似的数学验证 ====================

def mathematical_validation():
    """
    数学验证：Padé近似的连分式表示

    展示递归算法如何构造有理函数
    """
    print("\n" + "="*70)
    print("第3部分：Padé近似的数学结构")
    print("="*70)

    # 简单例子：用5个点构造Padé近似
    z_in = np.array([1j, 2j, 3j, 4j, 5j])
    f_in = 1.0 / (z_in + 10.0)  # 简单的解析函数

    # 构建g矩阵（手动展示递归过程）
    n = len(z_in)
    g = np.zeros((n, n), dtype=complex)
    g[0, :] = f_in

    print("Vidberg-Serene递归过程：")
    print("-" * 70)
    print("g矩阵的构建：")

    for i in range(1, n):
        for j in range(i, n):
            numerator = g[i-1, i-1] - g[i-1, j]
            denominator = (z_in[j] - z_in[i-1]) * g[i-1, j]
            g[i, j] = numerator / denominator

            if i < 3 and j < 3:  # 只打印前几个元素
                print(f"  g[{i},{j}] = {g[i,j]:.6f}")

    # 测试延拓
    z_test = 1.5 + 0.5j
    f_exact = 1.0 / (z_test + 10.0)
    f_pade = vidberg_serene_pade(z_in, f_in, z_test)

    print("\n延拓测试：")
    print(f"  测试点: z = {z_test}")
    print(f"  精确值: f(z) = {f_exact:.6f}")
    print(f"  Padé值: P(z) = {f_pade:.6f}")
    print(f"  相对误差: {abs(f_pade - f_exact)/abs(f_exact):.2e}")


# ==================== 主程序 ====================

def main():
    print("=" * 70)
    print("Padé解析延拓完整演示")
    print("=" * 70)

    # 第1部分：用已知函数验证
    test_pade_with_known_function()

    # 第2部分：Eliashberg实际应用
    example_eliashberg_continuation()

    # 第3部分：数学验证
    mathematical_validation()

    print("\n" + "=" * 70)
    print("总结:")
    print("-" * 70)
    print("1. Padé近似利用解析函数的唯一性从离散点重构函数")
    print("2. Vidberg-Serene算法通过递归构造有理函数")
    print("3. 对于光滑函数，Padé延拓精度很高（相对误差<1%）")
    print("4. Eliashberg求解器通过Padé延拓获得实验可观测量")
    print("5. 延拓质量依赖于Matsubara点的数量和分布")
    print("=" * 70)


if __name__ == "__main__":
    main()
