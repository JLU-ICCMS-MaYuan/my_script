#!/usr/bin/env python3
"""
将PROCAR文件转换为PBAND_La.dat格式
按照用户要求：遍历每个k点提取每个band的权重，将数据首尾相接
"""

import re
import sys
import numpy as np

def read_poscar(filename="POSCAR"):
    """读取POSCAR文件获取晶格参数"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 读取缩放因子
    scale = float(lines[1].strip())
    
    # 读取晶格向量
    lattice = []
    for i in range(2, 5):
        vec = [float(x) for x in lines[i].strip().split()]
        lattice.append(vec)
    
    lattice = np.array(lattice) * scale
    
    # 计算倒格子
    # 倒格子 = 2π * (正格子矩阵)^(-T)
    reciprocal_lattice = 2 * np.pi * np.linalg.inv(lattice).T
    
    return lattice, reciprocal_lattice

def parse_procar(filename):
    """解析PROCAR文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 从第二行获取k点数和band数
    header_line = lines[1].strip()
    match = re.search(r'# of k-points:\s*(\d+)\s*# of bands:\s*(\d+)', header_line)
    if match:
        nkpts = int(match.group(1))
        nbands = int(match.group(2))
    else:
        raise ValueError("无法解析k点数和band数")
    
    print(f"K点数: {nkpts}, Band数: {nbands}")
    
    # 存储数据的字典
    data = {}
    current_kpoint = None
    current_band = None
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 识别k点行
        if line.startswith('k-point'):
            match = re.search(r'k-point\s+(\d+)\s*:\s*([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)', line)
            if match:
                kpoint_idx = int(match.group(1)) - 1  # 转为0索引
                kx, ky, kz = float(match.group(2)), float(match.group(3)), float(match.group(4))
                current_kpoint = kpoint_idx
                if current_kpoint not in data:
                    data[current_kpoint] = {'kpoint': (kx, ky, kz), 'bands': {}}
        
        # 识别band行
        elif line.startswith('band'):
            match = re.search(r'band\s+(\d+)\s*#\s*energy\s*([\d\.\-]+)', line)
            if match:
                band_idx = int(match.group(1)) - 1  # 转为0索引
                energy = float(match.group(2))
                current_band = band_idx
                
                # 读取轨道权重数据：跳过空行和ion标题行
                i += 1  # 跳过空行
                if i < len(lines) and lines[i].strip() == '':
                    i += 1  # 跳过空行
                if i < len(lines) and lines[i].strip().startswith('ion'):
                    i += 1  # 跳过ion标题行
                
                if i < len(lines):
                    weight_line = lines[i].strip()
                    if weight_line.split()[0] == '1':  # ion 1的数据
                        weights = weight_line.split()[1:]  # 跳过离子序号
                        # 包含所有轨道：s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2, fy3x2, fxyz, fyz2, fz3, fxz2, fzx2, fx3, tot
                        orbital_weights = [float(w) for w in weights]  # 包含所有轨道权重
                        
                        if current_kpoint is not None:
                            data[current_kpoint]['bands'][current_band] = {
                                'energy': energy,
                                'weights': orbital_weights
                            }
        
        i += 1
    
    return data, nkpts, nbands

def calculate_kpath(kpoints, nkpts, reciprocal_lattice):
    """计算k路径距离，使用倒格子将分数坐标转换为绝对坐标"""
    kpath = [0.0]
    total_distance = 0.0
    
    for i in range(1, nkpts):
        prev_k_frac = np.array(kpoints[i-1]['kpoint'])
        curr_k_frac = np.array(kpoints[i]['kpoint'])
        
        # 将分数坐标转换为倒空间中的绝对坐标
        prev_k_cart = np.dot(prev_k_frac, reciprocal_lattice)
        curr_k_cart = np.dot(curr_k_frac, reciprocal_lattice)
        
        # 计算绝对距离
        distance = np.linalg.norm(curr_k_cart - prev_k_cart)
        total_distance += distance
        kpath.append(total_distance)
    
    return kpath

def write_pband_la(data, nkpts, nbands, output_filename, reciprocal_lattice):
    """写入PBAND_La.dat文件"""
    
    # 计算k路径
    kpoints = [data[i] for i in range(nkpts)]
    kpath = calculate_kpath(kpoints, nkpts, reciprocal_lattice)
    
    with open(output_filename, 'w') as f:
        # 写入头部
        f.write("#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot\n")
        f.write(f"# NKPTS & NBANDS: {nkpts} {nbands}\n")
        
        # 按照用户要求：遍历每个band，将所有k点的数据首尾相接
        for band_idx in range(nbands):
            f.write(f"# Band-Index    {band_idx + 1}\n")
            
            # 遍历所有k点，提取当前band的数据
            for kpoint_idx in range(nkpts):
                if kpoint_idx in data and band_idx in data[kpoint_idx]['bands']:
                    band_data = data[kpoint_idx]['bands'][band_idx]
                    energy = band_data['energy']
                    weights = band_data['weights']
                    
                    # 格式化输出：K-Path, Energy, 各轨道权重
                    line = f"{kpath[kpoint_idx]:10.5f}  {energy:12.6f}"
                    for w in weights:
                        line += f"  {w:5.3f}"
                    f.write(line + "\n")

def main():
    input_file = "PROCAR"
    output_file = f"PBAND_{sys.argv[1]}.dat"
    poscar_file = "POSCAR"
    
    print("开始读取POSCAR文件...")
    _, reciprocal_lattice = read_poscar(poscar_file)
    print("倒格子矩阵:")
    print(reciprocal_lattice)
    
    print("开始解析PROCAR文件...")
    data, nkpts, nbands = parse_procar(input_file)
    
    print(f"成功解析 {len(data)} 个k点")
    print("开始写入PBAND_La.dat文件...")
    
    write_pband_la(data, nkpts, nbands, output_file, reciprocal_lattice)
    
    print(f"转换完成！输出文件：{output_file}")

if __name__ == "__main__":
    main()