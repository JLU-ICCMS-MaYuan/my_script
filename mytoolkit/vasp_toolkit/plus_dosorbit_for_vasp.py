#!/usr/bin/env python3
import numpy as np
import pandas as pd
from ase.io import read
import argparse
from scipy import interpolate

# 命令行参数解析
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process PDOS data for multiple elements.")
    parser.add_argument('-e', '--elements', type=str, nargs='+',
                        help="List of elements with their orbitals, e.g., 'La:s,p,d,f' or 'H:s,p'.")
    parser.add_argument('-pv', '--per-volume', action='store_true',
                        help="Flag to indicate if the PDOS should be divided by the volume.")
    return parser.parse_args()

# 读取单个元素的 spdf 投影数据
def plus_spdf(e, orbitals_to_read):
    df = pd.read_table(f"PDOS_{e}.dat", sep=r'\s+')
    s = p = d = f = None

    if 's' in orbitals_to_read:
        s = df['s']
    if 'p' in orbitals_to_read:
        p = df['py'] + df['pz'] + df['px']
    if 'd' in orbitals_to_read:
        d = df['dxy'] + df['dyz'] + df['dz2'] + df['dxz'] + df['dx2']
    if 'f' in orbitals_to_read:
        f = (df['fy3x2'] + df['fxyz'] + df['fyz2'] + df['fz3'] +
             df['fxz2'] + df['fzx2'] + df['fx3'])
    tot = df['tot']
    return s, p, d, f, tot

# 插值获取 E=0 处的值
def interpolate_and_predict(energy, pdos):
    try:
        tck = interpolate.make_interp_spline(x=energy, y=pdos, k=1)
        return tck(0)
    except Exception:
        idx = np.abs(energy).argmin()
        return pdos.iloc[idx]

# 统一处理每个轨道
def process_orbital(orb, label, element, vol, per_volume, energy, result_df, pdos_at_ef):
    if orb is not None:
        column_name = f'{element}_{label}'
        if per_volume:
            result_df[column_name] = orb / vol
        else:
            result_df[column_name] = orb
        pdos_at_ef[column_name] = [interpolate_and_predict(energy, result_df[column_name])]

def main():
    args = parse_arguments()

    # 解析元素与轨道映射
    orbitals_dict = {}
    for item in args.elements:
        try:
            element, orbitals = item.split(':')
            orbitals_dict[element] = orbitals.split(',')
        except ValueError:
            raise ValueError(f"Invalid format for element input: {item}. Expected format like 'La:s,p,d'.")

    # 获取体积
    atoms = read('POSCAR')
    vol = atoms.get_volume()
    print(f'Volume of cell: {vol:.3f} Å³')

    # 读入参考 energy 列
    sample_element = list(orbitals_dict.keys())[0]
    df_sample = pd.read_table(f"PDOS_{sample_element}.dat", sep=r'\s+')
    energy = df_sample[df_sample.columns[0]]  # 自动适配 #Energy 或 Energy 等列名
    result_df = pd.DataFrame({'Energy(eV)': energy})
    pdos_at_ef = {'Energy(eV)': [0]}

    # 遍历元素处理
    for e, orbitals in orbitals_dict.items():
        s, p, d, f, tot = plus_spdf(e, orbitals)
        process_orbital(s, 's', e, vol, args.per_volume, energy, result_df, pdos_at_ef)
        process_orbital(p, 'p', e, vol, args.per_volume, energy, result_df, pdos_at_ef)
        process_orbital(d, 'd', e, vol, args.per_volume, energy, result_df, pdos_at_ef)
        process_orbital(f, 'f', e, vol, args.per_volume, energy, result_df, pdos_at_ef)

        colname = f'{e}'
        if args.per_volume:
            result_df[colname] = tot / vol
        else:
            result_df[colname] = tot
        pdos_at_ef[colname] = [interpolate_and_predict(energy, result_df[colname])]

    # 保存 Ef 插值值
    pdos_at_ef_df = pd.DataFrame(pdos_at_ef)
    pdos_at_ef_df.to_csv('PDOS_atEf.csv', index=False)

    # TDOS 插值
    tdos_df = pd.read_table("TDOS.dat", sep=r'\s+')
    energy_col = tdos_df.columns[0]
    energy_tdos = tdos_df[energy_col]
    tdos = tdos_df.iloc[:, 1]
    if args.per_volume:
        tdos_at_fermi = interpolate_and_predict(energy_tdos, tdos) / vol
        result_df['TDOS'] = tdos / vol
    else:
        tdos_at_fermi = interpolate_and_predict(energy_tdos, tdos)
        result_df['TDOS'] = tdos
    print(f"TDOS at Fermi level: {tdos_at_fermi:.6f}")

    # 保存总表
    output_file = 'spdf_orbit_per_volume.csv' if args.per_volume else 'spdf_orbit.csv'
    result_df.to_csv(output_file, index=False)
    print(f"Result saved to {output_file}")

    # 打印能量最接近 0 附近的数据
    idx_zero = result_df['Energy(eV)'].abs().idxmin()
    start = max(0, idx_zero - 3)
    end = min(len(result_df), idx_zero + 3)
    print("Around Fermi level (±3 rows):")
    print(result_df.iloc[start:end])

    # 计算 PDOS 总和 vs TDOS 差异
    columns_to_sum = [col for col in pdos_at_ef_df.columns if col not in ['Energy(eV)', 'TDOS']]
    pdos_sum_at_fermi = pdos_at_ef_df[columns_to_sum].sum(axis=1).item()
    diff = tdos_at_fermi - pdos_sum_at_fermi
    print(f"Sum of PDOS at Fermi level: {pdos_sum_at_fermi:.6f}")
    print(f"Difference (TDOS - sum of PDOS): {diff:.6f}")

if __name__ == "__main__":
    main()

