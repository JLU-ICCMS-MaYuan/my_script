#!/usr/bin/env python3
import numpy as np
import pandas as pd
from ase.io import read
import argparse
from scipy import interpolate

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process PDOS data for multiple elements")
    parser.add_argument('-e', '--elements', type=str, nargs='+', 
                        help="List of elements with their orbitals, e.g., 'La:s,p,d,f' or 'H:s,p'.")
    parser.add_argument('-pv', '--per-volume', action='store_true', default=True,
                        help="Flag to indicate if the PDOS should be divided by the volume.")
    return parser.parse_args()

def plus_spdf(e, orbitals_to_read):
    t_s_orb, t_p_orb, t_d_orb, t_f_orb, tot = [0, 0, 0, 0, 0]
    df = pd.read_table(f"PDOS_{e}.dat", sep='\s+')
    
    # 动态判断读取哪些轨道
    if 's' in orbitals_to_read:
        s_orb = df['s']
        t_s_orb += s_orb
    else:
        s_orb = None

    if 'p' in orbitals_to_read:
        p_orb = df['py'] + df['pz'] + df['px']
        t_p_orb += p_orb
    else:
        p_orb = None

    if 'd' in orbitals_to_read:
        d_orb = df['dxy'] + df['dyz'] + df['dz2'] + df['dxz'] + df['dx2']
        t_d_orb += d_orb
    else:
        d_orb = None

    if 'f' in orbitals_to_read:
        f_orb = df['fy3x2'] + df['fxyz'] + df['fyz2'] + df['fz3'] + df['fxz2'] + df['fzx2'] + df['fx3']
        t_f_orb += f_orb
    else:
        f_orb = None
    
    tot = df['tot']
    return t_s_orb, t_p_orb, t_d_orb, t_f_orb, tot

# 拟合并返回能量为0时的预测值，使用插值函数
def interpolate_and_predict(energy, pdos):
    tck = interpolate.make_interp_spline(x=energy, y=pdos, k=1)  # k=1 for linear spline
    return tck(0)  # 计算能量为0时的值

def main():
    # 解析命令行参数
    args = parse_arguments()

    # 创建一个字典，存储元素与轨道类型的映射
    orbitals_dict = {}
    for element_info in args.elements:
        element, orbitals = element_info.split(':')
        orbitals_dict[element] = orbitals.split(',')

    # 读取POSCAR文件并计算体积
    atoms = read('POSCAR')
    vol = atoms.get_volume()
    print('vol=', vol)

    eles = list(orbitals_dict.keys())  # 从字典中获取元素列表

    # 选择一个元素来获取 'Energy' 列
    sample_element = eles[0]
    sample_df = pd.read_table(f"PDOS_{sample_element}.dat", sep='\s+')
    energy = sample_df['#Energy']
    
    # 创建一个空 DataFrame 用于存储结果，并添加 'Energy' 列
    result_df = pd.DataFrame({'Energy(eV)': energy})

    # 用于存储能量为0时的PDOS值，初始化时使用 'Energy(eV)' 作为行索引
    pdos_at_ef = {'Energy(eV)': [0]}  # 初始化为能量列
    
    # 遍历元素
    for e in eles:
        orbitals_to_read = orbitals_dict[e]  # 获取该元素的轨道类型
        s_orb, p_orb, d_orb, f_orb, tot = plus_spdf(e, orbitals_to_read)
        
        # 将每个元素的 spdf 轨道和 tot 轨道添加到结果 DataFrame 中
        if s_orb is not None and args.per_volume:
            result_df[f'{e}_s'] = s_orb / vol
            pdos_at_ef[f'{e}_s'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_s'])]
        else:
            result_df[f'{e}_s'] = s_orb
            pdos_at_ef[f'{e}_s'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_s'])]
            
        if p_orb is not None and args.per_volume:
            result_df[f'{e}_p'] = p_orb / vol
            pdos_at_ef[f'{e}_p'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_p'])]
        else:
            result_df[f'{e}_p'] = p_orb
            pdos_at_ef[f'{e}_p'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_p'])]
            
        if d_orb is not None and args.per_volume:
            result_df[f'{e}_d'] = d_orb / vol
            pdos_at_ef[f'{e}_d'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_d'])]
        else:
            result_df[f'{e}_d'] = d_orb
            pdos_at_ef[f'{e}_d'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_d'])]

        if f_orb is not None and args.per_volume:
            result_df[f'{e}_f'] = f_orb / vol
            pdos_at_ef[f'{e}_f'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_f'])]
        else:
            result_df[f'{e}_f'] = f_orb
            pdos_at_ef[f'{e}_f'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}_f'])]
        
        if args.per_volume:    
            result_df[f'{e}'] = tot / vol
            pdos_at_ef[f'{e}'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}'])]
        else:
            result_df[f'{e}'] = tot
            pdos_at_ef[f'{e}'] = [interpolate_and_predict(result_df['Energy(eV)'], result_df[f'{e}'])]


    # 将能量为0时的PDOS值写入文件
    pdos_at_ef_df = pd.DataFrame(pdos_at_ef)
    pdos_at_ef_df.to_csv('PDOS_atEf.csv', index=False)

    # 获取费米能级处的 TDOS 值
    tdos_df = pd.read_table("TDOS.dat", sep='\s+')
    energy_tdos = tdos_df['#Energy']
    tdos = tdos_df['TDOS']
    # 使用插值函数插值计算能量为0时的 TDOS
    if args.per_volume:
        tdos_at_fermi = interpolate_and_predict(energy_tdos, tdos) / vol
        # 将 TDOS 数据添加到 result_df 最后一列
        result_df['TDOS'] = tdos/vol
    else:
        tdos_at_fermi = interpolate_and_predict(energy_tdos, tdos)
        # 将 TDOS 数据添加到 result_df 最后一列
        result_df['TDOS'] = tdos
    print(f"TDOS at Fermi level: {tdos_at_fermi}")
    
    # 保存结果到 spdf_orbit.csv
    result_df.to_csv('spdf_orbit.csv', index=False)
    print("Result saved to spdf_orbit.csv")
    nearest_to_zero_index = result_df['Energy(eV)'].abs().idxmin()
    # 获取最接近0的行的上下5行索引范围
    start_index = max(0, nearest_to_zero_index - 3)
    end_index = min(len(result_df), nearest_to_zero_index + 3)
    # 输出选取的行
    print("上下5行数据:")
    print(result_df.iloc[start_index:end_index])
    
    # 计算所有元素的 PDOS 值之和
    total_pdos_at_fermi = sum([pdos_at_ef_df[ele] for ele in eles])
    print(f"Sum of PDOS at Fermi level: {total_pdos_at_fermi.item()}")

    # 比较 TDOS 和 PDOS 之和
    diff = tdos_at_fermi - total_pdos_at_fermi
    print(f"Difference between TDOS and sum of PDOS: {diff.item()}")

if __name__ == "__main__":
    main()