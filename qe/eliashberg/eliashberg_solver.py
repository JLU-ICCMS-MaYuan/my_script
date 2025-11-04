#!/usr/bin/env python3
"""
Python版各向同性Eliashberg方程求解器 (从Fortran代码eliashberg.f90完整移植)

本脚本完全重现了Fortran版本的工作流程：

1. 从INPUT文件读取库仑赝势参数mu*和温度步数ntemp
2. 从ALPHA2F.OUT文件读取Eliashberg谱函数α²F(ω)
3. 计算McMillan/Allen-Dynes参数 (λ, ω_log, ω_rms, T_c)
4. 对每个温度在Matsubara轴上求解耦合的Eliashberg方程组：
   - Z(iωₙ) = 1 + (πT/ωₙ)∑ₘ (ωₘ/√(ωₘ²+Δₘ²)) λ(iωₙ-iωₘ)
   - Z(iωₙ)Δₙ = πT∑ₘ (Δₘ/√(ωₘ²+Δₘ²)) [λ(iωₙ-iωₘ) - μ*]
5. 使用Padé解析延拓将间隙函数和重整化函数延拓到实轴
6. 输出与Fortran版本完全一致的结果文件

数值细节(积分方案、混合参数、收敛容差)与参考实现完全匹配，
因此可以直接使用现有的输入文件。

核心物理：
- 基于Migdal-Eliashberg理论求解超导能隙
- 考虑电子-声子耦合和库仑排斥
- 适用于各向同性超导体的临界温度和能隙计算
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np


# 物理常数 (原子单位制Hartree)
PI = np.pi                    # π常数
TWO_PI = 2.0 * np.pi         # 2π，用于频率变换
KBOLTZ = 3.166815343e-6      # 玻尔兹曼常数 (Hartree/Kelvin)

# 迭代控制参数 (与Fortran代码完全相同)
MAX_WF = 40000               # 最大允许Matsubara频率数
MAX_IT = 1000                # 最大迭代次数
MIXING_BETA = 0.5            # 混合参数β (新旧解的线性组合系数)
CONV_TOL = 1.0e-12           # 收敛容差


@dataclass
class TemperatureSolution:
    """存储单个温度下Eliashberg方程的完整解"""
    temperature: float          # 温度 (Kelvin)
    matsubara_freq: np.ndarray  # Matsubara频率 ωₙ = πT(2n+1), 形状: (2*nwf+1,)
    gap_matsubara: np.ndarray   # Matsubara轴上的能隙函数 Δ(iωₙ)
    z_matsubara: np.ndarray     # Matsubara轴上的重整化函数 Z(iωₙ)
    gap_real_axis: np.ndarray   # 实轴上的能隙函数 Δ(ω) (复数)
    z_real_axis: np.ndarray     # 实轴上的重整化函数 Z(ω) (复数)
    real_freq: np.ndarray       # 实轴频率网格 (大小nout)
    iterations: int             # 达到收敛所需的迭代次数
    converged: bool             # 是否收敛标志
    nwf: int                    # 声子部分使用的Matsubara频率数
    nwf_coulomb: int            # 库仑部分使用的Matsubara频率数


def read_input(input_path: Path) -> Tuple[float, int]:
    """读取INPUT文件

    格式: mu* ntemp
    其中 mu* 是库仑赝势参数, ntemp 是温度步数
    这与Fortran版本的输入格式完全一致
    """
    line = input_path.read_text(encoding="ascii").splitlines()[0]
    tokens = line.split()
    if len(tokens) < 2:
        raise ValueError("INPUT 必须提供 mu* 与 ntemp 两个数值")
    mustar = float(tokens[0])    # 库仑赝势参数 μ*
    ntemp = int(tokens[1])       # 温度步数
    return mustar, ntemp


def read_alpha2f(alpha2f_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """读取Eliashberg谱函数α²F(ω)

    文件格式: 两列数据
    第一列: 频率 ω (原子单位)
    第二列: Eliashberg函数 α²F(ω)

    α²F(ω)包含了电子-声子耦合的全部信息，是求解超导能隙的核心输入
    """
    data = np.loadtxt(alpha2f_path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("ALPHA2F.OUT 需至少包含两列：频率 与 α²F(ω)")
    w = data[:, 0].astype(float)     # 频率网格
    a2f = data[:, 1].astype(float)   # Eliashberg谱函数
    idx = np.argsort(w)              # 按频率排序
    return w[idx], a2f[idx]


def trapezoid_integral(x: np.ndarray, y: np.ndarray) -> float:
    """梯形积分法

    用于计算∫ y(x) dx，与Fortran版本使用相同的数值积分方法
    """
    return float(np.trapz(y, x))


def compute_mcmillan_parameters(
    w: np.ndarray, a2f: np.ndarray, mustar: float
) -> Tuple[float, float, float, float, float, float]:
    """计算McMillan-Allen-Dynes参数

    基于Eliashberg谱函数计算关键的超导参数：
    1. λ: 电子-声子耦合常数 = 2∫ α²F(ω)/ω dω
    2. ω_log: 对数平均声子频率 = exp((2/λ)∫ α²F(ω)ln(ω)/ω dω)
    3. ω_rms: 均方根声子频率 = √((2/λ)∫ α²F(ω)ω dω)
    4. T_c: McMillan-Allen-Dynes超导临界温度公式

    这些参数与Fortran版本mcmillan.f90中的计算完全一致
    """
    # 避免除零，对极小频率做安全处理
    mask = w > 1.0e-12
    safe_w = np.where(mask, w, 1.0)

    # 计算电子-声子耦合常数 λ = 2∫ α²F(ω)/ω dω
    integrand_lambda = np.where(mask, a2f / safe_w, 0.0)
    lam = 2.0 * trapezoid_integral(w, integrand_lambda)

    # 计算对数平均频率 ω_log = exp((2/λ)∫ α²F(ω)ln(ω)/ω dω)
    integrand_log = np.where(mask, a2f * np.log(safe_w) / safe_w, 0.0)
    wlog = np.exp((2.0 / lam) * trapezoid_integral(w, integrand_log))

    # 计算均方根频率 ω_rms = √((2/λ)∫ α²F(ω)ω dω)
    integrand_rms = a2f * w
    wrms = np.sqrt(abs((2.0 / lam) * trapezoid_integral(w, integrand_rms)))

    # McMillan-Allen-Dynes超导临界温度公式
    # T_c = (ω_log/1.2k_B) * exp(-1.04(1+λ)/(λ-μ*-0.62λμ*)) * f1 * f2
    exponent = (-1.04 * (1.0 + lam)) / (lam - mustar - 0.62 * lam * mustar)
    tc = (wlog / (1.2 * KBOLTZ)) * np.exp(exponent)

    # Allen-Dynes修正因子 (Fig1.jpg形式)
    l1 = 2.46 * (1.0 + 3.8 * mustar)                    # 强耦合修正参数
    l2 = 1.82 * (1.0 + 6.3 * mustar) * (wrms / wlog)    # 声子谱形状修正参数
    f1 = (1.0 + (lam / l1) ** 1.5) ** (1.0 / 3.0)       # 强耦合修正因子
    f2 = 1.0 + (wrms / wlog - 1.0) * lam**2 / (lam**2 + l2**2)  # 谱形状修正因子
    tc *= f1 * f2

    return lam, wlog, wrms, tc, f1, f2


def compute_allen_dynes_comparison(
    w: np.ndarray, a2f: np.ndarray, mustar: float
) -> Tuple[dict, dict]:
    """计算两种Allen-Dynes修正形式的对比

    比较Fig1.jpg和Fig2.jpg中f2因子的两种不同表达式：

    Fig1形式 (Allen & Dynes 1975):
    f2 = 1 + λ²(ω̄₂/ωlog - 1) / [λ² + (1.82(1+6.3μ*)(ω̄₂/ωlog))²]

    Fig2形式 (简化形式):
    f1f2 = [1 + (λ/2.46(1+3.8μ*))^(3/2)]^(1/3) × [1 - λ²(1-ω₂/ωlog)/(λ² + 3.312(1+6.3μ*)²)]

    返回两种形式的详细参数对比
    """
    # 首先计算基本参数
    mask = w > 1.0e-12
    safe_w = np.where(mask, w, 1.0)

    # 电子-声子耦合常数
    integrand_lambda = np.where(mask, a2f / safe_w, 0.0)
    lam = 2.0 * trapezoid_integral(w, integrand_lambda)

    # 对数平均频率
    integrand_log = np.where(mask, a2f * np.log(safe_w) / safe_w, 0.0)
    wlog = np.exp((2.0 / lam) * trapezoid_integral(w, integrand_log))

    # 均方根频率
    integrand_rms = a2f * w
    wrms = np.sqrt(abs((2.0 / lam) * trapezoid_integral(w, integrand_rms)))

    # 基础临界温度 (McMillan指数部分)
    exponent = (-1.04 * (1.0 + lam)) / (lam - mustar - 0.62 * lam * mustar)
    tc_base = (wlog / (1.2 * KBOLTZ)) * np.exp(exponent)

    # Fig1形式 (Allen & Dynes 1975原始公式)
    fig1_results = {}
    fig1_results['form'] = 'Fig1 (Allen & Dynes 1975)'
    fig1_results['lambda'] = lam
    fig1_results['omega_log'] = wlog * 315.77e3  # 转换为cm⁻¹
    fig1_results['omega_rms'] = wrms * 315.77e3  # 转换为cm⁻¹
    fig1_results['omega_ratio'] = wrms / wlog

    # f1因子 (两种形式相同)
    l1 = 2.46 * (1.0 + 3.8 * mustar)
    f1 = (1.0 + (lam / l1) ** 1.5) ** (1.0 / 3.0)
    fig1_results['f1'] = f1
    fig1_results['l1'] = l1

    # f2因子 (Fig1形式)
    omega_ratio = wrms / wlog
    l2_fig1 = 1.82 * (1.0 + 6.3 * mustar) * omega_ratio
    f2_fig1 = 1.0 + (omega_ratio - 1.0) * lam**2 / (lam**2 + l2_fig1**2)
    fig1_results['f2'] = f2_fig1
    fig1_results['l2'] = l2_fig1
    fig1_results['tc_kelvin'] = tc_base * f1 * f2_fig1

    # Fig2形式 (简化形式)
    fig2_results = {}
    fig2_results['form'] = 'Fig2 (Simplified form)'
    fig2_results['lambda'] = lam
    fig2_results['omega_log'] = wlog * 315.77e3  # 转换为cm⁻¹
    fig2_results['omega_rms'] = wrms * 315.77e3  # 转换为cm⁻¹
    fig2_results['omega_ratio'] = wrms / wlog

    # f1因子 (相同)
    fig2_results['f1'] = f1
    fig2_results['l1'] = l1

    # f2因子 (Fig2形式) - 注意这里是组合的f1f2因子
    l2_fig2_squared = 3.312 * (1.0 + 6.3 * mustar)**2
    f2_correction = 1.0 - lam**2 * (1.0 - omega_ratio) / (lam**2 + l2_fig2_squared)
    f1f2_combined = f1 * f2_correction
    fig2_results['f2_correction'] = f2_correction
    fig2_results['f1f2_combined'] = f1f2_combined
    fig2_results['l2_squared'] = l2_fig2_squared
    fig2_results['tc_kelvin'] = tc_base * f1f2_combined

    return fig1_results, fig2_results


def build_lambda_kernel(
    w: np.ndarray, a2f: np.ndarray, t0: float, nmax: int
) -> np.ndarray:
    """构建电子-声子耦合核 Λ_m

    计算Eliashberg方程中的耦合核函数:
    Λ_m = 2∫ [ωα²F(ω) / (ω² + (2t_0 m)²)] dω

    这对应Fortran代码中的l(−m)的计算，其中:
    - t_0 = πk_B T 是Matsubara温度
    - m ∈ [-2nmax, 2nmax] 是频率索引差

    参数:
        w: 声子频率网格
        a2f: Eliashberg谱函数 α²F(ω)
        t0: Matsubara温度 πk_B T
        nmax: 最大Matsubara频率索引

    返回:
        形状为(4*nmax+1,)的数组，存储所有可能的Λ_m值
    """
    m_vals = np.arange(-2 * nmax, 2 * nmax + 1, dtype=int)
    kernel = np.empty_like(m_vals, dtype=float)
    w_sq = w**2
    numerator = w * a2f

    # 对每个频率差 m 计算耦合核
    for idx, m in enumerate(m_vals):
        denom = w_sq + (2.0 * t0 * m) ** 2  # 分母: ω² + (2t_0 m)²
        # 避免除零，对分母为0的情况返回0
        integrand = np.divide(numerator, denom, out=np.zeros_like(w_sq), where=denom > 0)
        kernel[idx] = 2.0 * trapezoid_integral(w, integrand)
    return kernel


def pade_approximation(
    zin: np.ndarray, uin: np.ndarray, zout: np.ndarray
) -> np.ndarray:
    """使用Vidberg & Serene递归算法进行Padé解析延拓

    将Matsubara轴(iωₙ)上的函数延拓到实轴(ω)上。
    这是计算实频率谱函数的关键步骤。

    算法来源: H. J. Vidberg and J. W. Serene,
    J. Low Temp. Phys. 29, 179 (1977)

    参数:
        zin: 输入复数频率点 (Matsubara频率 iωₙ)
        uin: 输入函数值 (能隙或Z函数在Matsubara轴上的值)
        zout: 输出复数频率点 (实轴频率)

    返回:
        延拓到实轴上的函数值 (复数)

    这个实现与Fortran版本pade.f90完全一致
    """
    nin = zin.size
    nout = zout.size
    if nin == 0 or nout == 0:
        return np.zeros(nout, dtype=complex)

    # 按照Vidberg-Serene方法构建g矩阵
    g = np.zeros((nin, nin), dtype=complex)
    g[0, :] = uin  # 初始化第一行

    # 递归构建g矩阵的其余元素
    for i in range(1, nin):
        for j in range(i, nin):
            denom = (zin[j] - zin[i - 1]) * g[i - 1, j]
            g[i, j] = (g[i - 1, i - 1] - g[i - 1, j]) / denom

    # 对每个输出点计算Padé近似
    uout = np.zeros(nout, dtype=complex)
    for i, z in enumerate(zout):
        # 初始化递归参数
        a0 = 0.0 + 0.0j
        a1 = g[0, 0]
        b0 = 1.0 + 0.0j
        b1 = 1.0 + 0.0j

        # 递归计算Padé近似的分子和分母
        for j in range(1, nin):
            zt1 = (z - zin[j - 1]) * g[j, j]
            a0, a1 = a1, a1 + zt1 * a0  # 更新分子项
            b0, b1 = b1, b1 + zt1 * b0  # 更新分母项

            # 防止数值溢出，必要时进行缩放
            if max(abs(a1.real), abs(a1.imag)) > 1.0e100:
                scale = 1.0 / abs(a1)
                a0 *= scale
                a1 *= scale
                b0 *= scale
                b1 *= scale
            if max(abs(b1.real), abs(b1.imag)) > 1.0e100:
                scale = 1.0 / abs(b1)
                a0 *= scale
                a1 *= scale
                b0 *= scale
                b1 *= scale

        # 计算最终结果 a1/b1
        if abs(b1.real) + abs(b1.imag) > 0.0:
            uout[i] = a1 / b1
    return uout


def fold_matsubara_vectors(values: np.ndarray) -> np.ndarray:
    """构建对称的Matsubara频率数组

    Eliashberg方程求解时只计算正频率部分(ωₙ ≥ 0)，
    然后利用对称性扩展到负频率部分:
    - 对于能隙: Δ(-ωₙ) = Δ(ωₙ) (偶函数)
    - 对于Z函数: Z(-ωₙ) = Z(ωₙ) (偶函数)

    这对应Fortran代码中的对称性处理逻辑

    参数:
        values: 正频率部分的值 (大小: nwf+1)

    返回:
        完整的对称数组 (大小: 2*nwf+1)
    """
    nmax = values.size - 1  # nwf
    folded = np.empty(2 * nmax + 1, dtype=float)

    # 对每个索引 n ∈ [-nmax, nmax] 赋值
    for idx, n in enumerate(range(-nmax, nmax + 1)):
        if n >= 0:
            m = n           # 正频率: 直接使用
        else:
            m = -n - 1      # 负频率: 利用对称性
        folded[idx] = values[m]
    return folded


def solve_eliashberg(
    mustar: float,
    ntemp: int,
    w: np.ndarray,
    a2f: np.ndarray,
    lam: float,
    wlog: float,
    wrms: float,
    tc: float,
) -> List[TemperatureSolution]:
    """求解各向同性Eliashberg方程组

    这是整个程序的核心函数，在给定的温度范围内求解耦合的Eliashberg方程:

    Z(iωₙ) = 1 + (πT/ωₙ)∑ₘ (ωₘ/√(ωₘ²+Δₘ²)) [Λ(ωₙ-ωₘ) - Λ(ωₙ+ωₘ)]
    Z(iωₙ)Δₙ = πT∑ₘ (Δₘ/√(ωₘ²+Δₘ²)) [Λ(ωₙ-ωₘ) + Λ(ωₙ+ωₘ)] - 2μ*∑ₘ (Δₘ Zₘ/√(ωₘ²+Δₘ²))

    其中:
    - ωₙ = πT(2n+1) 是Matsubara频率
    - Δₙ 是超导能隙函数
    - Zₙ 是准粒子重整化函数
    - Λ(ω) 是电子-声子耦合核
    - μ* 是库仑赝势参数

    计算流程:
    1. 设置温度范围和频率截断
    2. 对每个温度点进行迭代求解
    3. 使用混合算法加速收敛
    4. 进行Padé解析延拓到实轴

    返回每个温度下的完整解，包括能隙和Z函数在Matsubara轴和实轴上的值
    """
    # 设置计算参数 (与Fortran版本完全一致)
    wfmax = 20.0 * wrms      # Matsubara频率截断: 20ω_rms
    tmin = tc / 6.0          # 最低温度: T_c/6
    tmax = 3.0 * tc          # 最高温度: 3T_c
    dtemp = (tmax - tmin) / float(ntemp)  # 温度步长

    # 分配内存和初始化
    nwf_alloc = int(round(wfmax / (TWO_PI * KBOLTZ * dtemp)))  # 预估最大频率数
    nwf_alloc = min(max(nwf_alloc, 1), MAX_WF)               # 限制在合理范围内
    d0 = np.full(nwf_alloc + 1, 1.0e-4, dtype=float)        # 能隙初始值 (1e-4)
    z0 = np.ones(nwf_alloc + 1, dtype=float)                # Z函数初始值 (1.0)

    solutions: List[TemperatureSolution] = []
    # 主循环: 遍历所有温度点
    for itemp in range(1, ntemp + 1):
        temp = itemp * dtemp + tmin      # 当前温度
        t0 = PI * KBOLTZ * temp          # Matsubara温度参数 πk_B T

        # 计算当前温度下的频率数
        nwf = int(round(wfmax / (2.0 * t0)))     # 声子部分频率数
        nwf = min(max(nwf, 1), nwf_alloc)      # 限制在分配范围内
        nwfcl = int(round(20.0 * wrms / (2.0 * t0)))  # 库仑部分频率数
        nwfcl = max(1, min(nwfcl, nwf))         # 保证不超过声子部分

        # 生成Matsubara频率 ωₙ = πT(2n+1)
        wf_pos = t0 * (2 * np.arange(0, nwf + 1, dtype=float) + 1.0)
        # 计算电子-声子耦合核
        lambda_kernel = build_lambda_kernel(w, a2f, t0, nwf)
        offset = 2 * nwf  # 数组索引偏移量

        # 复制初始值作为当前温度的起始点
        d_slice = d0[: nwf + 1].copy()  # 能隙函数 Δ(iωₙ)
        z_slice = z0[: nwf + 1].copy()  # Z函数 Z(iωₙ)

        # 自洽迭代求解Eliashberg方程
        converged = False
        iterations = 0
        for it in range(1, MAX_IT + 1):
            iterations = it

            # 计算准粒子谱函数的分母: r = |Zω| = √((Zω)² + (ZΔ)²)
            r = np.sqrt((wf_pos**2 + d_slice**2) * z_slice**2)
            r = np.where(r == 0.0, 1.0e-30, r)  # 避免除零

            # 求解Z函数方程: Z(iωₙ) = 1 + (πT/ωₙ)∑ₘ (ωₘ/rₘ) [Λ(ωₙ-ωₘ) - Λ(ωₙ+ωₘ)]
            z_new = np.empty_like(z_slice)
            for n in range(nwf + 1):
                total = 0.0
                for m in range(nwf):
                    # 获取耦合核元素
                    lm = lambda_kernel[(n - m) + offset]      # Λ(ωₙ-ωₘ)
                    lp = lambda_kernel[(n + m + 1) + offset]  # Λ(ωₙ+ωₘ)
                    # 累加求和项
                    total += (lm - lp) * z_slice[m] * wf_pos[m] / r[m]
                z_new[n] = t0 * total / wf_pos[n]  # 乘以 πT/ωₙ
            z_new += 1.0       # 加上常数项
            z_slice = z_new

            # 计算库仑修正项: 2μ*∑ₘ (Δₘ Zₘ/rₘ)
            # 只在低频率(nwfcl)部分考虑库仑效应
            dmu = 2.0 * mustar * np.sum(
                d_slice[: nwfcl + 1] * z_slice[: nwfcl + 1] / r[: nwfcl + 1]
            )

            # 求解能隙方程: Z(iωₙ)Δₙ = πT∑ₘ (Δₘ/rₘ) [Λ(ωₙ-ωₘ) + Λ(ωₙ+ωₘ)] - dmu
            d_new = np.empty_like(d_slice)
            for n in range(nwf + 1):
                total = 0.0
                for m in range(nwf):
                    # 获取耦合核元素
                    lm = lambda_kernel[(n - m) + offset]      # Λ(ωₙ-ωₘ)
                    lp = lambda_kernel[(n + m + 1) + offset]  # Λ(ωₙ+ωₘ)
                    # 累加求和项 (注意这里是加号)
                    total += (lm + lp) * d_slice[m] * z_slice[m] / r[m]
                # 除以Z函数得到能隙
                d_new[n] = t0 * (total - dmu) / z_slice[n]

            # 线性混合加速收敛: d_new = β*d_new + (1-β)*d_old
            mixed = MIXING_BETA * d_new + (1.0 - MIXING_BETA) * d_slice
            # 计算收敛判据: 平均绝对差值
            diff = np.sum(np.abs(d_slice - mixed)) / max(1, 2 * nwf)
            d_slice = mixed

            # 检查是否收敛
            if diff <= CONV_TOL:
                converged = True
                break

        # 保存收敛结果作为下一个温度的初始值
        d0[: nwf + 1] = d_slice
        z0[: nwf + 1] = z_slice

        # 构建完整的对称Matsubara频率数组 [-nwf, nwf]
        freq_full = t0 * (2 * np.arange(-nwf, nwf + 1, dtype=float) + 1.0)
        gap_full = fold_matsubara_vectors(d_slice)   # 对称扩展能隙函数
        z_full = fold_matsubara_vectors(z_slice)     # 对称扩展Z函数

        # 准备Padé解析延拓到实轴
        denom = max(len(w) - 1, 1)
        dw = (w[-1] - w[0]) / float(denom)           # 频率步长
        nout = 4 * len(w)                            # 输出点数
        real_freq = 2.0 * np.arange(nout, dtype=float) * dw  # 实轴频率网格

        # 执行Padé解析延拓
        zin = 1j * wf_pos.astype(complex)            # Matsubara频率 iωₙ
        # 将能隙和Z函数从Matsubara轴延拓到实轴
        gap_real = pade_approximation(zin, d_slice.astype(complex), real_freq + 0j)
        z_real = pade_approximation(zin, z_slice.astype(complex), real_freq + 0j)

        solutions.append(
            TemperatureSolution(
                temperature=temp,
                matsubara_freq=freq_full,
                gap_matsubara=gap_full,
                z_matsubara=z_full,
                gap_real_axis=gap_real,
                z_real_axis=z_real,
                real_freq=real_freq,
                iterations=iterations,
                converged=converged,
                nwf=nwf,
                nwf_coulomb=nwfcl,
            )
        )
    return solutions


def write_allen_dynes_comparison(
    output_dir: Path,
    fig1_results: dict,
    fig2_results: dict
) -> None:
    """输出Allen-Dynes两种形式的对比结果

    生成AllenDynes_Tcs.dat文件，包含两种f2修正形式的详细对比
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ad_path = output_dir / "AllenDynes_Tcs.dat"

    with ad_path.open("w", encoding="utf-8") as fh:
        fh.write("+" + "="*60 + "+\n")
        fh.write("|" + " Allen-Dynes f2因子两种形式对比分析 ".center(60) + "|\n")
        fh.write("+" + "="*60 + "+\n\n")

        # 基本参数
        fh.write("基本参数:\n")
        fh.write("-" * 40 + "\n")
        fh.write(f"λ (电子-声子耦合常数)    : {fig1_results['lambda']:12.6f}\n")
        fh.write(f"ω_log (对数平均频率)     : {fig1_results['omega_log']:12.2f} cm⁻¹\n")
        fh.write(f"ω_rms (均方根频率)       : {fig1_results['omega_rms']:12.2f} cm⁻¹\n")
        fh.write(f"ω_rms/ω_log (频率比)     : {fig1_results['omega_ratio']:12.6f}\n")
        fh.write(f"f1 (强耦合修正因子)      : {fig1_results['f1']:12.6f}\n\n")

        # Fig1形式详细结果
        fh.write("Fig1形式 (Allen & Dynes 1975原始公式):\n")
        fh.write("-" * 50 + "\n")
        fh.write("公式: f2 = 1 + λ²(ω̄₂/ωlog - 1) / [λ² + (1.82(1+6.3μ*)(ω̄₂/ωlog))²]\n\n")
        fh.write(f"l1 = 2.46(1+3.8μ*)       : {fig1_results['l1']:12.6f}\n")
        fh.write(f"l2 = 1.82(1+6.3μ*)(ω̄₂/ωlog): {fig1_results['l2']:12.6f}\n")
        fh.write(f"f2 (谱形状修正因子)      : {fig1_results['f2']:12.6f}\n")
        fh.write(f"T_c (临界温度)           : {fig1_results['tc_kelvin']:12.4f} K\n\n")

        # Fig2形式详细结果
        fh.write("Fig2形式 (简化形式):\n")
        fh.write("-" * 40 + "\n")
        fh.write("公式: f1f2 = [1 + (λ/l1)^(3/2)]^(1/3) × [1 - λ²(1-ω₂/ωlog)/(λ² + 3.312(1+6.3μ*)²)]\n\n")
        fh.write(f"l2² = 3.312(1+6.3μ*)²    : {fig2_results['l2_squared']:12.6f}\n")
        fh.write(f"f2修正项                 : {fig2_results['f2_correction']:12.6f}\n")
        fh.write(f"f1×f2 (组合修正因子)     : {fig2_results['f1f2_combined']:12.6f}\n")
        fh.write(f"T_c (临界温度)           : {fig2_results['tc_kelvin']:12.4f} K\n\n")

        # 差异分析
        tc_diff = fig1_results['tc_kelvin'] - fig2_results['tc_kelvin']
        tc_rel_diff = abs(tc_diff) / fig1_results['tc_kelvin'] * 100
        f2_diff = fig1_results['f2'] - fig2_results['f2_correction']

        fh.write("差异分析:\n")
        fh.write("-" * 30 + "\n")
        fh.write(f"T_c差异                  : {tc_diff:+12.4f} K\n")
        fh.write(f"T_c相对差异              : {tc_rel_diff:12.3f} %\n")
        fh.write(f"f2修正因子差异           : {f2_diff:+12.6f}\n\n")

        # 公式转换说明
        fh.write("公式关系说明:\n")
        fh.write("-" * 40 + "\n")
        fh.write("两种形式在物理上等价，差异主要来自：\n")
        fh.write("1. Fig1: 保留完整的(ω̄₂/ωlog)²项在分母中\n")
        fh.write("2. Fig2: 将(ω̄₂/ωlog)²≈1近似，简化分母\n")
        fh.write("3. Fig2: 将1-ω₂/ωlog写成-(ω₂/ωlog-1)形式\n")
        fh.write(f"4. 当ω̄₂/ωlog={fig1_results['omega_ratio']:.3f}时，近似误差为{tc_rel_diff:.3f}%\n\n")

        # 推荐使用
        if tc_rel_diff < 1.0:
            fh.write("推荐: 两种形式差异很小(<1%)，可任选其一\n")
        elif fig1_results['omega_ratio'] > 1.2 or fig1_results['omega_ratio'] < 0.8:
            fh.write("推荐: ω̄₂/ωlog偏离1较多，建议使用Fig1完整形式\n")
        else:
            fh.write("推荐: 可使用Fig2简化形式，但Fig1更精确\n")

    print(f"Allen-Dynes对比结果已写入: {ad_path}")


def write_output_files(
    output_dir: Path,
    mustar: float,
    ntemp: int,
    w: np.ndarray,
    lam: float,
    wlog: float,
    wrms: float,
    tc: float,
    f1: float,
    f2: float,
    solutions: List[TemperatureSolution],
) -> None:
    """输出结果文件

    生成与Fortran版本完全一致的输出文件:
    1. ELIASHBERG.OUT: 主要计算信息和参数总结
    2. ELIASHBERG_IA.OUT: Matsubara轴上的能隙和Z函数
    3. ELIASHBERG_GAP_T.OUT: 能隙随温度的变化
    4. ELIASHBERG_GAP_RA.OUT: 实轴上的能隙函数
    5. ELIASHBERG_Z_RA.OUT: 实轴上的Z函数

    所有文件格式与Fortran版本保持一致，便于对比和后续分析
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    nout = 4 * len(w)
    wfmax = 20.0 * wrms
    tmin = tc / 6.0
    tmax = 3.0 * tc

    with (output_dir / "ELIASHBERG.OUT").open("w", encoding="ascii") as fh:
        fh.write("+------------------------------+\n")
        fh.write("|     Eliashberg equations     |\n")
        fh.write("+------------------------------+\n\n")
        fh.write(f"lambda : {lam:18.10E}\n")
        fh.write(f"omega_log : {wlog:18.10E}\n")
        fh.write(f"omega_rms : {wrms:18.10E}\n")
        fh.write(f"f1 : {f1:18.10E}\n")
        fh.write(f"f2 : {f2:18.10E}\n")
        fh.write(f"mu* : {mustar:18.10E}\n\n")
        fh.write(f"Temperature range : {tmin:18.10E} {tmax:18.10E}\n")
        fh.write(f"Number of temperature steps : {ntemp:6d}\n")
        fh.write(f"Number of output frequencies : {nout:8d}\n")
        fh.write("Fermionic Matsubara frequency cut-off\n")
        fh.write(f" phonons : {wfmax:18.10E}\n")
        fh.write(f" Coulomb : {wrms:18.10E}\n")
        for sol in solutions:
            fh.write("\n")
            fh.write(f"Temperature (kelvin) : {sol.temperature:18.10E}\n")
            fh.write("Number of Matsubara frequencies\n")
            fh.write(f" phonons : {sol.nwf:8d}\n")
            fh.write(f" Coulomb : {sol.nwf_coulomb:8d}\n")
            if sol.converged:
                fh.write(
                    f"Eliashberg equations converged in {sol.iterations:6d} iterations\n"
                )
            else:
                fh.write("Failed to converge: possibly close to T_c\n")

    ia_path = output_dir / "ELIASHBERG_IA.OUT"
    with ia_path.open("w", encoding="ascii") as fh:
        for sol in solutions:
            for freq, gap, zval in zip(
                sol.matsubara_freq, sol.gap_matsubara, sol.z_matsubara
            ):
                fh.write(f"{freq:18.10E} {gap:18.10E} {zval:18.10E}\n")
            fh.write("\n")

    gap_t_path = output_dir / "ELIASHBERG_GAP_T.OUT"
    with gap_t_path.open("w", encoding="ascii") as fh:
        for sol in solutions:
            fh.write(
                f"{sol.temperature:18.10E} {sol.gap_matsubara[sol.nwf]:18.10E} "
                f"{sol.z_matsubara[sol.nwf]:18.10E}\n"
            )

    gap_ra_path = output_dir / "ELIASHBERG_GAP_RA.OUT"
    with gap_ra_path.open("w", encoding="ascii") as fh:
        for sol in solutions:
            for freq, val in zip(sol.real_freq, sol.gap_real_axis):
                fh.write(
                    f"{freq:18.10E} {val.real:18.10E} {val.imag:18.10E}\n"
                )
            fh.write("\n")

    z_ra_path = output_dir / "ELIASHBERG_Z_RA.OUT"
    with z_ra_path.open("w", encoding="ascii") as fh:
        for sol in solutions:
            for freq, val in zip(sol.real_freq, sol.z_real_axis):
                fh.write(
                    f"{freq:18.10E} {val.real:18.10E} {val.imag:18.10E}\n"
                )
            fh.write("\n")

    print()
    print("Info(eliashberg):")
    print(f" calculation information written to {ia_path.parent/'ELIASHBERG.OUT'}")
    print(
        " gap and Z functions on the imaginary axis written to "
        f"{ia_path}"
    )
    print(
        " gap vs. temperature written to "
        f"{gap_t_path}"
    )
    print(
        " gap function on the real axis written to "
        f"{gap_ra_path}"
    )
    print(
        " Z function on the real axis written to "
        f"{z_ra_path}"
    )


def main() -> None:
    """主函数: 各向同性Eliashberg方程求解器

    执行完整的计算流程:
    1. 解析命令行参数
    2. 读取输入文件 (INPUT, ALPHA2F.OUT)
    3. 计算McMillan-Allen-Dynes参数
    4. 求解耦合的Eliashberg方程
    5. 输出结果文件

    使用方法:
        python3 eliashberg_solver.py [--input INPUT] [--alpha2f ALPHA2F.OUT] [--output-dir .]
    """
    parser = argparse.ArgumentParser(
        description="求解各向同性Eliashberg方程 (Python版本)"
    )
    parser.add_argument(
        "--input", default="INPUT", type=str, help="包含 mu* 与 ntemp 的文件"
    )
    parser.add_argument(
        "--alpha2f", default="ALPHA2F.OUT", type=str, help="Eliashberg 谱函数文件"
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        type=str,
        help="输出目录（默认当前目录）",
    )
    args = parser.parse_args()

    # 准备文件路径
    input_path = Path(args.input)
    alpha_path = Path(args.alpha2f)
    output_dir = Path(args.output_dir)

    # 执行主计算流程
    mustar, ntemp = read_input(input_path)                        # 读取基本参数
    w, a2f = read_alpha2f(alpha_path)                            # 读取谱函数
    lam, wlog, wrms, tc, f1, f2 = compute_mcmillan_parameters(w, a2f, mustar)  # 计算特征参数

    # 计算Allen-Dynes两种形式的对比
    fig1_results, fig2_results = compute_allen_dynes_comparison(w, a2f, mustar)
    write_allen_dynes_comparison(output_dir, fig1_results, fig2_results)

    solutions = solve_eliashberg(mustar, ntemp, w, a2f, lam, wlog, wrms, tc)    # 求解方程
    write_output_files(                                         # 输出结果
        output_dir, mustar, ntemp, w, lam, wlog, wrms, tc, f1, f2, solutions
    )


if __name__ == "__main__":
    main()
