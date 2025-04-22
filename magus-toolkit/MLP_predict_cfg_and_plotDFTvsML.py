#!/usr/bin/env python3
import subprocess
import argparse

import numpy as np
from ase.atoms import Atoms
from ase.units import GPa
from collections import defaultdict

import matplotlib.pyplot as plt

def dump_cfg(frames, filename, symbol_to_type, mode='w'):
    with open(filename, mode) as f:
        for atoms in frames:
            ret = ''
            ret += 'BEGIN_CFG\n'
            ret += 'Size\n{}\n'.format(len(atoms))
            try:
                cell = atoms.get_cell()[:]
                ret += 'Supercell\n{} {} {}\n{} {} {}\n{} {} {}\n'.format(*cell[0], *cell[1], *cell[2])
            except:
                pass
            cartes = atoms.positions
            has_forces = False
            has_constraints = False
            try:
                atoms.info['forces'] = atoms.get_forces()
            except:
                pass
            fields = ['id', 'type', 'cartes_x', 'cartes_y', 'cartes_z']
            if 'forces' in atoms.info:
                fields.extend(['fx', 'fy', 'fz'])
                forces = atoms.info['forces']
                has_forces = True

            flags = np.array(['T']*len(atoms))
            if atoms.constraints:
                atoms.info.append('mvable')
                for constr in atoms.constraints:
                    flags[constr.index] = 'F'
                has_constraints = True

            ret += 'AtomData: ' + ' '.join(fields) + '\n'
            for i, atom in enumerate(atoms):
                atom_info = '{} {} {} {} {} '.format(i + 1, symbol_to_type[atom.symbol], *cartes[i])
                if has_forces:
                    atom_info += '{} {} {} '.format(*forces[i])
                if has_constraints:
                    atom_info += flags[i]
                ret += atom_info + '\n'
            try:
                atoms.info['energy'] = atoms.get_potential_energy()
            except:
                pass
            if 'energy' in atoms.info:
                ret += 'Energy\n{}\n'.format(atoms.info['energy'])
            try:
                atoms.info['stress'] = atoms.get_stress()
            except:
                pass
            if 'stress' in atoms.info:
                stress = atoms.info['stress'] * atoms.get_volume() * -1.
                ret += 'PlusStress: xx yy zz yz xz xy\n{} {} {} {} {} {}\n'.format(*stress)
            if 'identification' in atoms.info:
                ret += 'Feature identification {}\n'.format(atoms.info['identification'])
            ret += 'END_CFG\n'
            f.write(ret)


#TODO 
# different cell
# pbc
def load_cfg(filename, type_to_symbol):
    frames = []
    with open(filename) as f:
        line = 'chongchongchong!'
        while line:
            line = f.readline()
            if 'BEGIN_CFG' in line:
                cell = np.zeros((3, 3))

            if 'Size' in line:
                line = f.readline()
                natoms = int(line.split()[0])
                positions = np.zeros((natoms, 3))
                forces = np.zeros((natoms, 3))
                energies = np.zeros(natoms)
                symbols = ['X'] * natoms
                constraints = ['T'] * natoms

            if 'Supercell' in line: 
                for i in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[i, j] = float(line.split()[j])

            if 'AtomData' in line:
                d = defaultdict(int)
                for (i, x) in enumerate(line.split()[1:]):
                    d[x] = i

                for _ in range(natoms):
                    line = f.readline()
                    fields = line.split()
                    i = int(fields[d['id']]) - 1
                    symbols[i] = type_to_symbol[int(fields[d['type']])]
                    positions[i] = [float(fields[d[attr]]) for attr in ['cartes_x', 'cartes_y' ,'cartes_z']]
                    forces[i] = [float(fields[d[attr]]) for attr in ['fx', 'fy' ,'fz']]
                    energies[i] = float(fields[d['site_en']])

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                if d['fx'] != 0:
                    atoms.info['forces'] = forces
                if d['site_en'] != 0:
                    atoms.info['energies'] = energies

            if 'Energy' in line and 'Weight' not in line:
                line = f.readline()
                atoms.info['energy'] = float(line.split()[0])

            if 'PlusStress' in line:
                line = f.readline()
                plus_stress = np.array(list(map(float, line.split())))
                atoms.info['PlusStress'] = plus_stress
                atoms.info['stress'] = -plus_stress / atoms.get_volume()
                atoms.info['pstress'] = atoms.info['stress'] / GPa

            if 'END_CFG' in line:
                frames.append(atoms)

            if 'identification' in line:
                atoms.info['identification'] = int(line.split()[2])

    return frames


# 提取能量并存储
def extract_energies(DFTframes, MLPframes):
    DFTenergies = []
    MLPenergies = []
    for DFT, MLP in zip(DFTframes, MLPframes):
        # energy1 = DFT.info.get('energy', None)
        # energy2 = MLP.info.get('energy', None)
        assert len(DFT) == len(MLP)
        DFTenergies.append(DFT.info['energy']/len(DFT))
        MLPenergies.append(MLP.info['energy']/len(MLP))
    return np.array(DFTenergies), np.array(MLPenergies)


# 提取力并存储
def extract_forces(DFTframes, MLPframes):
    DFTforces = []
    MLPforces = []
    for DFT, MLP in zip(DFTframes, MLPframes):
        # print(DFT.symbols, MLP.symbols)
        # print(DFT.get_scaled_positions(), MLP.get_scaled_positions())
        # print(DFT.info.get('forces', None))
        # print(MLP.info.get('forces', None))
        DFTforces.append(DFT.info['forces'])
        MLPforces.append(MLP.info['forces'])
    return DFTforces, MLPforces


# 绘制能量对比图
def plot_energy_comparison(DFTenergies, MLPenergies):
    plt.figure(figsize=(8, 6))
    plt.scatter(DFTenergies, MLPenergies, color='blue')
    plt.xlabel('DFT energies (eV)')
    plt.ylabel('MLP energies (eV)')
    plt.title('Energy Comparison')
    
    # 添加对角线（淡灰色）
    max_energy = max(max(DFTenergies), max(MLPenergies))
    min_energy = min(min(DFTenergies), min(MLPenergies))
    plt.plot([min_energy, max_energy], [min_energy, max_energy], color='red', linestyle='--', alpha=0.5)
    
    plt.grid(True)
    plt.savefig('energy_comparison.png')  # 保存图像
    plt.close()


# # 绘制力对比图
# def plot_force_comparison(DFTforces, MLPforces):
#     fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
#     directions = ['x', 'y', 'z']

#     for k, ax in enumerate(axes):  # 遍历 x, y, z 分量
#         DFT_flatten = []
#         MLP_flatten = []

#         for DFTforce, MLPforce in zip(DFTforces, MLPforces):
#             DFT_flatten.extend(DFTforce[:, k])  # 提取第 k 列 (x, y, z 分量)
#             MLP_flatten.extend(MLPforce[:, k])

#         # 绘制散点图
#         ax.scatter(DFT_flatten, MLP_flatten, color='blue', alpha=0.5)
#         ax.set_xlabel(f'DFT Force {directions[k]} (eV/Å)')
#         ax.set_ylabel(f'MLP Force {directions[k]} (eV/Å)')
#         ax.set_title(f'Force Comparison ({directions[k]})')
#         ax.grid(True)

#         # 添加对角线
#         xrange = np.linspace(min(DFT_flatten), max(DFT_flatten))
#         ax.plot(xrange, xrange, color='red', linestyle='--', alpha=0.7)

#     plt.tight_layout()
#     plt.savefig('force_comparison.png')  # 保存图像
#     plt.close()                
         
# 绘制力对比图
def plot_force_comparison(DFTforces, MLPforces):
    directions = ['x', 'y', 'z']
    colors = ['blue', 'green', 'orange']
    
    plt.figure(figsize=(8, 8))
    diag_line = []
    for k in range(3):  # 遍历 x, y, z 分量
        DFT_flatten = []
        MLP_flatten = []
        
        for DFTforce, MLPforce in zip(DFTforces, MLPforces):
            DFT_flatten.extend(DFTforce[:, k])  # 提取第 k 列 (x, y, z 分量)
            MLP_flatten.extend(MLPforce[:, k])
    
        # 绘制散点图
        plt.scatter(DFT_flatten, MLP_flatten, color=colors[k], alpha=0.5, label=f'Force {directions[k]}')
        diag_line.extend([min(DFT_flatten), max(DFT_flatten), min(MLP_flatten), max(MLP_flatten)])
        
    # 添加对角线
    force_range =  np.linspace(min(diag_line), max(diag_line))
    plt.plot(force_range, force_range, color='red', linestyle='--', alpha=0.7)

    # 图形设置
    plt.xlabel('DFT Force (eV/Å)')
    plt.ylabel('MLP Force (eV/Å)')
    plt.title('Force Comparison (x, y, z)')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('force_comparison.png')  # 保存图像
    plt.close()
    
    
def find_energy_deviation(DFTenergies:np.ndarray, MLPenergies:np.ndarray, threshold=1.0):
    """
    使用 NumPy 数据广播找到能量偏差超过指定阈值的结构索引号。

    参数：
    - DFTenergies: list 或 numpy array, DFT计算得到的能量。
    - MLPenergies: list 或 numpy array, MLP计算得到的能量。
    - threshold: float, 偏差阈值, 默认是1.0 eV。

    返回：
    - indices: numpy array, 偏差小于等于阈值的结构索引号。
    """
    deviations = np.abs(MLPenergies - DFTenergies)
    indices = np.where(deviations <= threshold)[0]  # 筛选满足条件的索引
    return indices


def filter_force_deviation(DFTforces:list[np.ndarray], MLPforces:np.ndarray, threshold=200):
    """
    根据力的偏差筛选结构。

    参数：
        - structures: list, 包含多个结构数据。
        - threshold: float, 筛选阈值, 单位 eV/Å。

    返回：
        - 筛选后的结构列表。
    """
    indices = []
    i = 0
    for DFTforce, MLPforce in zip(DFTforces, MLPforces):
        deviations = np.abs(MLPforce - DFTforce)
        if np.all(deviations < threshold):
            indices.append(i)
            # print(np.max(deviations))
        i+=1
    return indices


def save_filtered_cfg(DFTframes, indices, output_filename, symbol_to_type):
    """
    根据筛选条件保存符合偏差条件的结构到新的CFG文件。

    参数：
    - DFTframes: list of Atoms, 原始的DFT计算的结构数据。
    - indices: list 或 numpy array, 满足条件的结构索引号。
    - output_filename: str, 保存的CFG文件名。
    - symbol_to_type: dict, 元素符号到类型的映射关系。
    """
    filtered_frames = [DFTframes[i] for i in indices]
    dump_cfg(filtered_frames, output_filename, symbol_to_type)


def execute_command(command, output_file=None):
    """
    执行命令并将输出写入文件（可选）。
    - command: str, 要执行的命令。
    - output_file: str or None, 输出重定向的文件名。如果为 None, 则不保存输出。
    """
    with open(output_file, 'w') if output_file else subprocess.DEVNULL as f_out:
        process = subprocess.Popen(
            command, 
            shell=True, 
            stdout=f_out, 
            stderr=subprocess.STDOUT
        )
        process.wait()
        
               
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("read_cfg.py -dft train.cfg  -mlp calculated.cfg   -s Ce Sr H -t 500 10 -o train1.cfg")
    parser.add_argument('-dft', '--dftcfg',  default='train.cfg',   help='please input cfg file name')
    parser.add_argument('-mlp', '--mlpcfg', default='calculated.cfg',   help='please input cfg file name')
    parser.add_argument('-s',   '--symbols',   default='', nargs='+', help='name of elements')
    parser.add_argument('-t',   '--threshold', type=float, nargs='+', help='能量偏差的阈值, 单位为eV, 力偏差阈值, 单位eV/A')
    parser.add_argument('-o',   '--output', default='new_train.cfg', help='保存筛选结果的CFG文件名')
    parser.add_argument('-e',   '--execute', action="store_true", help="如果指定此选项, 则实际执行mlp命令, 否则仅打印命令")
    args = parser.parse_args()
    
    # 构造命令
    calc_efs_cmd = f"mlp calc-efs pot.mtp {args.dftcfg} {args.mlpcfg}"
    calc_errors_cmd = f"mlp calc-errors pot.mtp {args.dftcfg}"

    # 执行命令并将输出保存到文件
    calc_efs_output = "calc_efs.log"
    calc_errors_output = "calc_errors.log"
    
    # 打印并根据命令行参数控制是否执行 calc_efs_cmd
    print(calc_efs_cmd)
    print(f"输出保存到 {calc_efs_output}")
    # 打印并根据命令行参数控制是否执行 calc_errors_cmd
    print(calc_errors_cmd)
    print(f"输出保存到 {calc_errors_output}")
    if args.execute:
        execute_command(calc_efs_cmd, calc_efs_output)
        execute_command(calc_errors_cmd, calc_errors_output)

    
    type_to_symbol = {i: j for i, j in enumerate(args.symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]
    
    # 载入数据, 并判断dftcfg的结构数是否等于mlpcfg的结构数
    DFTframes = load_cfg(args.dftcfg, type_to_symbol)
    MLPframes = load_cfg(args.mlpcfg, type_to_symbol)
    assert len(DFTframes) == len(MLPframes), f"Length mismatch: len(DFTframes) = {len(DFTframes)}, len(MLPframes) = {len(MLPframes)}"

    # 提取能量和力数据
    DFTenergies, MLPenergies = extract_energies(DFTframes, MLPframes) # eV/atom
    DFTforces, MLPforces = extract_forces(DFTframes, MLPframes)
    
    # 绘制能量对比图
    plot_energy_comparison(DFTenergies, MLPenergies)
    
    # 绘制力对比图
    plot_force_comparison(DFTforces, MLPforces)

    # 筛选偏差小于等于阈值的结构
    indices1 = find_energy_deviation(DFTenergies, MLPenergies, threshold=args.threshold[0])
    indices2 = filter_force_deviation(DFTforces, MLPforces, threshold=args.threshold[1])
    print(f"筛选到 {len(indices1)} 个符合偏差小于等于 {args.threshold[0]} eV   的结构。")
    print(f"筛选到 {len(indices2)} 个符合偏差小于等于 {args.threshold[1]} eV/A 的结构。")
    indices = set(indices1) & set(indices2) 
    print(f"筛选到 {len(indices)} 个")
    # 保存筛选结果
    save_filtered_cfg(DFTframes, indices, args.output, symbol_to_type)
    print(f"筛选后的结构已保存到 {args.output}")