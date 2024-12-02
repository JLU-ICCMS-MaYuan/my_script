#!/usr/bin/env python3
import sys
import os
from pathlib import Path


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_eliashberg_Tc(eliashberg_x_path, idx):
    pass

def eliashberg_Tc(eliashberg_x_path, idx):
    mu_tc = {}
    for mu in [0.1, 0.13, 0.16]:
        with open('INPUT', 'w') as f:
            f.write("{:<10} {:<10}".format(mu, 1000))
        tc = get_eliashberg_Tc(eliashberg_x_path, idx)
        mu_tc[mu]=tc
    return mu_tc

if __name__ == "__main__":
    Ry2Thz = 3289.84195885094

    try:
        el_ph_sigma  = float(sys.argv[1]) # 起始展宽和展宽步长
    except:
        el_ph_sigma  = 0.005 # 起始展宽和展宽步长

    try:
        el_ph_nstep  = el_ph_sigma
    except:
        el_ph_nstep  = el_ph_sigma

    try:
        el_ph_nsigma = int(sys.argv[2]) # 展宽数
    except:
        el_ph_nsigma = 10 # 展宽数

    try:
        eliashberg_x_path = sys.argv[3]
    except:
        eliashberg_x_path = '~/code/my_script/qe/eliashberg/eliashberg.x'


    total_sigma  = np.linspace(el_ph_sigma, el_ph_sigma + (el_ph_nsigma - 1)*el_ph_nstep, num=el_ph_nsigma)
    print(total_sigma)

    integral_results = []; mu1_tc = []; mu2_tc = []; mu3_tc = []
    for idx in range(1, len(total_sigma)+1):
        print(f'sigma={total_sigma[idx-1]}')
        data = pd.read_csv('a2F.dos'+str(idx), delimiter='\s+', skiprows=5, header=None)
        # 删除最后一行
        data = data.iloc[:-1]
        data = data.apply(pd.to_numeric, errors='coerce')
        frequency = data.iloc[:,0]*Ry2Thz
        a2F_tot   = data.iloc[:,1]

        #print(frequency, a2F_tot)
        integrand = 2*a2F_tot/frequency
        integral = np.trapz(integrand, frequency)
        integral_results.append(integral)
        #print(f"Integral result for a2F.dos{idx}: {integral}")

        mu_tc = eliashberg_Tc(eliashberg_x_path, idx) # mu_tc = {0.1: tc1,  0.13: tc2,  0.16: tc3}
        mu1_tc.append(mu_tc[0.1])
        mu2_tc.append(mu_tc[0.13])
        mu3_tc.append(mu_tc[0.16])

    # 绘制图表
    # Plot integral_results with left y-axis
    fig, ax1 = plt.subplots()
    ax1.plot(total_sigma, integral_results, marker='o', linestyle='-', color='b', label='Integral Results')
    ax1.set_xlabel('Total Sigma')
    ax1.set_ylabel('Integral Results', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)  # 添加网格
    
    # Create a second y-axis
    ax2 = ax1.twinx()
    ax2.plot(total_sigma, mu1_tc, marker='o', linestyle='-', color='r', label='Mu1 = 0.1')
    ax2.plot(total_sigma, mu2_tc, marker='o', linestyle='-', color='g', label='Mu2 = 0.13')
    ax2.plot(total_sigma, mu3_tc, marker='o', linestyle='-', color='m', label='Mu3 = 0.16')
    ax2.set_ylabel('Mu TC', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.grid(True)  # 添加网格
    
    # Adding legends
    fig.tight_layout()
    fig.legend(loc='upper left', bbox_to_anchor=(0.1, 0.9))
    
    plt.show()
    plt.savefig('sigma-lambda-Tc.png')

    # 打印所有结果
    print("All integral results:", integral_results)

    # 创建 DataFrame
    df = pd.DataFrame({
        'Total Sigma': total_sigma,
        'Integral Results': integral_results,
        'Tc(mu=0.1)': mu1_tc,
        'Tc(mu=0.13)': mu2_tc,
        'Tc(mu=0.16)': mu3_tc,
    })
    df.to_csv('sigma_lambda_Tc.csv', index=False)

