#!/usr/bin/env python3

import sys
import pandas as pd

from ase.io import read

atoms = read('POSCAR')
vol = atoms.get_volume()
print('vol=', vol)

eles = sys.argv[1:]

df = pd.read_table("PDOS_H.dat", sep='\s+')
H_tdos = df['tot']
energy = df['#Energy']

df = pd.read_table("PDOS_Be.dat", sep='\s+')
Be_tdos = df['tot']

t_d_orb = 0
t_f_orb = 0
t_s_orb = 0
t_p_orb = 0
# 遍历重元素
for e in eles:
    df = pd.read_table(f"PDOS_{e}.dat", sep='\s+')
    s_orb = df['s'] 
    p_orb = df['py'] + df['pz'] + df['px']
    d_orb = df['dxy'] + df['dyz'] + df['dz2'] + df['dxz'] + df['dx2']
    f_orb = df['fy3x2'] + df['fxyz'] + df['fyz2'] + df['fz3'] + df['fxz2'] + df['fzx2'] + df['fx3'] 
    t_s_orb += s_orb
    t_p_orb += p_orb
    t_d_orb += d_orb
    t_f_orb += f_orb

# 创建新的DataFrame
new_df = pd.DataFrame({
    '#Energy': energy,
    '_s,'.join(eles)+'_s,'+'_p,'.join(eles)+'_p': (t_s_orb + t_p_orb)/vol,
    '_d,'.join(eles)+'_d': t_d_orb/vol,
    '_f,'.join(eles)+'_f': t_f_orb/vol,
    'Be': Be_tdos/vol,
    'H': H_tdos/vol,
    'TDOS': (H_tdos + Be_tdos + t_f_orb + t_d_orb + t_s_orb + t_p_orb)/vol,
})
nearest_to_zero_index = new_df['#Energy'].abs().idxmin()
# 获取最接近0的行的上下5行索引范围
start_index = max(0, nearest_to_zero_index - 3)
end_index = min(len(new_df), nearest_to_zero_index + 3)

# 输出选取的行
print("上下5行数据:")
print(new_df.iloc[start_index:end_index])

new_df.to_csv('final.csv',index=False)
