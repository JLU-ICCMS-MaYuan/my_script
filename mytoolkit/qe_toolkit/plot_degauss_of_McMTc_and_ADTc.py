#!/usr/bin/env python3
import sys
import os
from pathlib import Path


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_a2Fdos_data(gauss_idx):
    a2Fdos_path = Path("a2F.dos"+str(gauss_idx))
    if not a2Fdos_path.exists():
        print("\nNote: --------------------")
        print(f"    `a2F.dos+{str(gauss_idx)}` file doesn't exist!!!")
        sys.exit(1)
    # 这里简单说明一下a2F.dos*的文件结构
    # 频率             总a2F   剩下的列是N*3列，是N个原子，每个原子3个振动模式。
    # 特别注意第一列频率的单位是：里德伯Ry
    a2Fdos_data = os.popen(f"sed '1,5d' {a2Fdos_path} | sed '/lambda/d' ").readlines()
    a2Fdos_data = np.array([[float(x) for x in line.strip('\n').split()] for line in a2Fdos_data])

    return a2Fdos_data

def get_lambda_from_a2fdos_single_broadening(a2Fdos_data, gauss_idx):
    """
    输入: 
        gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
    """

    # 计算lambda
    frequency = a2Fdos_data[:,0]
    a2F = a2Fdos_data[:,1]
    lambda_value = np.trapz(2 * a2F / frequency, frequency)

    # 将 frequency 和 a2F 作为列来创建 DataFrame
    w_alpha2f = pd.DataFrame({
        'omegas(Rydberg)': frequency,
        'alpha2f': a2F
    })
    w_alpha2f.to_csv(
        "w_alpha2f_from_a2Fdos"+str(gauss_idx)+".csv",
        header=True,
        index=False,
    )
    return lambda_value

def get_wlog_from_a2fdos_single_broadening(a2Fdos_data, Lambda_bya2Fdos, gauss_idx):
    """
    输入: 
        gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
    """

    # 计算wlog
    Ry2K = 157887.51240116
    frequency = a2Fdos_data[:,0]
    a2F = a2Fdos_data[:,1]
    w_log = Ry2K * np.exp(2/Lambda_bya2Fdos * np.trapz(a2F / frequency * np.log(frequency), frequency))
    return w_log

def get_w2_from_a2fdos_single_broadening(a2Fdos_data, Lambda_bya2Fdos, gauss_idx):
    """
    输入: 
        gauss_idx: 是一个收敛的degauss的索引, 因为python是从0开始的, 所以一定要小心
        gauss: 对应索引的Gaussian值
    """
    # 计算w2
    frequency = a2Fdos_data[:,0]
    a2F = a2Fdos_data[:,1]
    w2  = np.sqrt(2/Lambda_bya2Fdos * np.trapz(frequency*a2F, frequency))
    return w2

def get_Tc_McM_Tc_AD(
            Lambda,
            wlog,
            w2):
    
    screen_constants = [0.1, 0.13, 0.16]
    Tcs_McM = {}; Tcs_AD = {}
    for screen_constant in screen_constants:
        #print(f"sed '1,5d' "+ "  a2F.dos"+str(idx) + "  | sed '/lambda/d' | awk '{print $1/2, $2}' > ALPHA2F.OUT")
        f1 = np.cbrt(1 + np.power(Lambda / (2.46 * (1 + 3.8 * screen_constant)), 1.5) )
        # f2 = 1 + (Lambda**2 * (1 - w2 / wlog)) / (Lambda**2 + 3.312 * ((1 + 6.3 * screen_constant) * w2 / wlog)**2)
        f2 = 1 - (Lambda**2 * (1-w2/wlog)) / (Lambda**2 + 3.312*(1+6.3*screen_constant)**2)
        Tc_McM = wlog/1.2 * np.exp( (-1.04*(1+Lambda)) / (Lambda-screen_constant*(1+0.62*Lambda)) )
        Tc_AD  = f1*f2*Tc_McM
        Tcs_McM[screen_constant] = Tc_McM
        Tcs_AD[screen_constant] = Tc_AD
    return Tcs_McM, Tcs_AD


if __name__ == "__main__":
    Ry2Thz = 3289.84195885094

    try:
        el_ph_sigma  = float(sys.argv[1]) # 起始展宽和展宽步长
    except:
        el_ph_sigma  = 0.005 # 起始展宽

    try:
        el_ph_nstep  = el_ph_sigma # 展开递增的步长
    except:
        el_ph_nstep  = el_ph_sigma

    try:
        el_ph_nsigma = int(sys.argv[2]) # 展宽数
    except:
        el_ph_nsigma = 10 # 展宽数


    total_sigma  = np.linspace(el_ph_sigma, el_ph_sigma + (el_ph_nsigma - 1)*el_ph_nstep, num=el_ph_nsigma)
    print(total_sigma)

    integral_results = []; tc_McM_mu1 = []; tc_McM_mu2 = []; tc_McM_mu3 = []; tc_AD_mu1 = []; tc_AD_mu2 = []; tc_AD_mu3 = []
    for idx in range(1, len(total_sigma)+1):
        print(f'sigma={total_sigma[idx-1]}')
        data = get_a2Fdos_data(idx)
        Lambda = get_lambda_from_a2fdos_single_broadening(data, idx)
        wlog = get_wlog_from_a2fdos_single_broadening(data, Lambda, idx)
        w2 = get_w2_from_a2fdos_single_broadening(data, Lambda, idx)
        tcs_McM, tcs_AD = get_Tc_McM_Tc_AD(Lambda, wlog, w2) # tcs_McM = {0.1: tc1,  0.13: tc2,  0.16: tc3}
        tc_McM_mu1.append(tcs_McM[0.1])
        tc_McM_mu2.append(tcs_McM[0.13])
        tc_McM_mu3.append(tcs_McM[0.16])
        tc_AD_mu1.append(tcs_AD[0.1])
        tc_AD_mu2.append(tcs_AD[0.13])
        tc_AD_mu3.append(tcs_AD[0.16])

    # 绘制图表
    plt.figure(figsize=(10, 6))
    plt.plot(total_sigma, tc_McM_mu1, label='Tc_McM (mu=0.1)',  marker='o')
    plt.plot(total_sigma, tc_McM_mu2, label='Tc_McM (mu=0.13)', marker='o')
    plt.plot(total_sigma, tc_McM_mu3, label='Tc_McM (mu=0.16)', marker='o')
    plt.xlabel('Lambda')
    plt.ylabel('Tc_McM (K)')
    plt.title('Tc_McM vs Lambda')
    plt.legend()
    plt.grid(True)
    plt.savefig('Tc_McM_vs_Lambda.png')  # 保存图像
    plt.show()

    # 2. 绘制 tc_AD_mu1, tc_AD_mu2, tc_AD_mu3 与 Lambda 的关系
    plt.figure(figsize=(10, 6))
    plt.plot(total_sigma, tc_AD_mu1, label='Tc_AD (mu=0.1)', marker='o')
    plt.plot(total_sigma, tc_AD_mu2, label='Tc_AD (mu=0.13)', marker='o')
    plt.plot(total_sigma, tc_AD_mu3, label='Tc_AD (mu=0.16)', marker='o')
    plt.xlabel('Lambda')
    plt.ylabel('Tc_AD (K)')
    plt.title('Tc_AD vs Lambda')
    plt.legend()
    plt.grid(True)
    plt.savefig('Tc_AD_vs_Lambda.png')  # 保存图像
    plt.show()

    # 打印所有结果
    print("All integral results:", integral_results)

    # 创建 DataFrame
    df = pd.DataFrame({
        'Total Sigma': total_sigma,
        'Integral Results': integral_results,
        'Tc(mu=0.1)': tc_McM_mu1,
        'Tc(mu=0.13)': tc_McM_mu2,
        'Tc(mu=0.16)': tc_McM_mu3,
    })
    df.to_csv('sigma_lambda_McM.csv', index=False)

    # 创建 DataFrame
    df = pd.DataFrame({
        'Total Sigma': total_sigma,
        'Integral Results': integral_results,
        'Tc(mu=0.1)': tc_AD_mu1,
        'Tc(mu=0.13)': tc_AD_mu2,
        'Tc(mu=0.16)': tc_AD_mu3,
    })
    df.to_csv('sigma_lambda_AD.csv', index=False)