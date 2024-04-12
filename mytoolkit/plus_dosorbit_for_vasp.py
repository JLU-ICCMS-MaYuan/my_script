#!/usr/bin/env python3

import sys
import pandas as pd

from ase.io import read


def plus_spdf(e):
    t_s_orb, t_p_orb, t_d_orb, t_f_orb, tot = [0,0,0,0,0]
    df = pd.read_table(f"PDOS_{e}.dat", sep='\s+')
    s_orb = df['s'] 
    p_orb = df['py'] + df['pz'] + df['px']
    d_orb = df['dxy'] + df['dyz'] + df['dz2'] + df['dxz'] + df['dx2']
    f_orb = df['fy3x2'] + df['fxyz'] + df['fyz2'] + df['fz3'] + df['fxz2'] + df['fzx2'] + df['fx3'] 
    t_s_orb += s_orb
    t_p_orb += p_orb
    t_d_orb += d_orb
    t_f_orb += f_orb
    tot = df['tot']
    return t_s_orb, t_p_orb, t_d_orb, t_f_orb, tot

atoms = read('POSCAR')
vol = atoms.get_volume()
print('vol=', vol)

eles = sys.argv[1:]

# 选择一个元素来获取 'Energy' 列
sample_element = eles[0]
sample_df = pd.read_table(f"PDOS_{sample_element}.dat", sep='\s+')
energy = sample_df['#Energy']
# 创建一个空 DataFrame 用于存储结果，并添加 'Energy' 列
result_df = pd.DataFrame({'#Energy': energy})


# 遍历元素
for e in eles:
    df = pd.read_table(f"PDOS_{e}.dat", sep='\s+')
    s_orb, p_orb, d_orb, f_orb, tot = plus_spdf(e)
    
    # 将每个元素的 spdf 轨道和 tot 轨道添加到结果 DataFrame 中
    result_df[f'{e}_s'] = s_orb/vol
    result_df[f'{e}_p'] = p_orb/vol
    result_df[f'{e}_d'] = d_orb/vol
    result_df[f'{e}_f'] = f_orb/vol
    result_df[f'{e}_tot'] = tot/vol

# 将结果写入文件
result_df.to_csv('spdf_tot_result.csv', index=False)

nearest_to_zero_index = result_df['#Energy'].abs().idxmin()
# 获取最接近0的行的上下5行索引范围
start_index = max(0, nearest_to_zero_index - 3)
end_index = min(len(result_df), nearest_to_zero_index + 3)

# 输出选取的行
print("上下5行数据:")
print(result_df.iloc[start_index:end_index])

result_df.to_csv('final.csv',index=False)
