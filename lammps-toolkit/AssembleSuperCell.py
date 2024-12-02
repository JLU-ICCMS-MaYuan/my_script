#!/usr/bin/env python3
from ase import Atoms
from ase.io import read, write
from ase.io.lammpsdata import write_lammps_data
from ase.build import make_supercell
from ase.calculators.lammps import Prism, convert

import numpy as np


def write_lammps_dump_text(
        filename: str,
        image: Atoms,
        specorder: list = None,
        write_image_flags: bool = False,
        reduce_cell: bool = False,
        units: str = 'metal',
        prismobj: Prism = None,
        ):
    symbols = image.get_chemical_symbols()
    if specorder is None:
        # This way it is assured that LAMMPS atom types are always
        # assigned predictably according to the alphabetic order
        species = sorted(set(symbols))
    else:
        # To index elements in the LAMMPS data file
        # (indices must correspond to order in the potential file)
        species = specorder
        
    with open(filename, 'w') as fd:
        fd.write("ITEM: TIMESTEP\n")
        fd.write("0\n")
        fd.write("ITEM: NUMBER OF ATOMS\n")
        fd.write(f"{len(image)}\n")
        fd.write("ITEM: BOX BOUNDS pp pp pp\n")
        if prismobj is None:
            prismobj = Prism(image.get_cell(), reduce_cell=reduce_cell)
        # Get cell parameters and convert from ASE units to LAMMPS units
        xhi, yhi, zhi, xy, xz, yz = convert(
            prismobj.get_lammps_prism(), 'distance', 'ASE', units)
        fd.write(f'0.0 {xhi:23.17g}\n')
        fd.write(f'0.0 {yhi:23.17g}\n')
        fd.write(f'0.0 {zhi:23.17g}\n')
        fd.write("ITEM: ATOMS id element xu yu zu\n")
        
        if write_image_flags:
            scaled_positions = image.get_scaled_positions(wrap=False)
            image_flags = np.floor(scaled_positions).astype(int)

        pos = prismobj.vector_to_lammps(
            image.get_positions(),
            wrap=write_image_flags,
        )
        
        pos = convert(pos, 'distance', 'ASE', tounits='metal')
        for i, r in enumerate(pos):
            s = species.index(symbols[i]) + 1
            element_name = species[s - 1]
            line = (
                # f'{i+1:>6} {s:>3}'
                    f'{i+1:>6} {element_name:>4}'
                f' {r[0]:23.17g} {r[1]:23.17g} {r[2]:23.17g}'
            )
            if write_image_flags:
                img = image_flags[i]
                line += f' {img[0]:6d} {img[1]:6d} {img[2]:6d}'
            line += '\n'
            fd.write(line)

def remove_hydrogen_and_convert(frame:Atoms, element_order):
    """
    读取 LAMMPS 轨迹文件，删除氢原子，按照指定顺序重新排列元素，写入 LAMMPS 格式文件。
    
    Args:
        input_file (str): 输入的 LAMMPS trajectory 文件名。
        output_file (str): 输出的 LAMMPS 格式文件名。
        element_order (list): 指定的元素排列顺序，例如 ['C', 'O', 'H', 'N']。
    """

    # 处理每一帧，删除氢原子并按照指定顺序排列元素
    symbols = frame.get_chemical_symbols()
    # 只保留非氢和非铍原子
    indices_to_keep = [i for i, symbol in enumerate(symbols) if symbol not in ['H', 'Be']]
    
    # 删除氢和铍原子
    new_frame = frame[indices_to_keep]
    
    # 按照指定的元素顺序排序
    element_to_index = {element: idx for idx, element in enumerate(element_order)}
    
    # 获取排序后的索引
    sorted_indices = sorted(range(len(new_frame)), key=lambda i: element_to_index[new_frame.get_chemical_symbols()[i]])
    
    # 根据指定顺序排列原子
    new_frame = new_frame[sorted_indices]

    return new_frame

def replace_element(structure, original_element, new_element):
    positions = structure.positions
    species = structure.get_chemical_symbols()
    # print("species", species); input()
    new_species = []
    for i, el in enumerate(species):
        if el == original_element:
            # 循环替换元素
            new_species.append(new_element)
        else:
            new_species.append(el)

    # print(new_species)
    structure.set_chemical_symbols(new_species)
    structure.positions = positions

    return structure

def stack_structures_along_c_axis(structures):
    """
    将多个结构沿着 c 轴方向叠加
    """
    all_positions = []
    all_symbols = []
    cell = structures[0].cell  # 假设所有结构的晶格相同
    c_length = cell[2, 2]  # 获取 c 轴长度

    for idx, structure in enumerate(structures):
        positions = structure.positions
        symbols = structure.get_chemical_symbols()

        # 将每个结构的原子位置沿着 c 轴方向平移
        offset = idx * c_length
        new_positions = positions + np.array([0, 0, offset])
        
        all_positions.append(new_positions)
        all_symbols.extend(symbols)

    # 合并所有的位置
    all_positions = np.vstack(all_positions)

    # 更新晶格：沿 c 轴扩展
    new_cell = structures[0].cell.copy()
    new_cell[2, 2] = c_length * len(structures)  # 更新 c 轴长度

    # 创建新的结构
    new_structure = Atoms(symbols=all_symbols, positions=all_positions, cell=new_cell, pbc=True)

    return new_structure


if __name__ == "__main__":
    
    # 1. 读取POSCAR文件
    poscar_file = 'La.vasp'
    structure = read(poscar_file)

    # 2. 扩胞 4x4x1
    cell = structure.cell
    # 定义扩胞矩阵
    supercell_matrix = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]])  # 4*4*1扩胞矩阵
    supercell = make_supercell(structure, supercell_matrix)

    # 3. 替换La为Y、Ce和Th
    old_elements = ['La']
    new_elements = ['Y', 'Ce', 'Th']
    allsupercell = [supercell]
    for old_ele in old_elements:
        for new_ele in new_elements:
            new_supercell = replace_element(supercell.copy(), old_ele, new_ele)
            #  new_supercell = replace_element(supercell, old_ele, new_ele) # 必须要用.copy()
            allsupercell.append(new_supercell)
    
    # 4. 沿c轴叠加所有结构
    stacked_structure = stack_structures_along_c_axis(allsupercell)
    
    # 5. 排序元素顺序
    # 获取所有的元素符号
    symbols = stacked_structure.get_chemical_symbols()
    element_order = ['La', 'Y', 'Ce', 'Th', 'Be', 'H']
    sorted_indices = sorted(range(len(stacked_structure)), key=lambda i: element_order.index(symbols[i]))
    
    # 6. 根据指定顺序排列原子
    stacked_structure = stacked_structure[sorted_indices]
    
    # 7. 输出新的结构
    write_lammps_data('stacked_structure.nso', stacked_structure, specorder=['La', 'Y', 'Ce', 'Th', 'Be', 'H'])
    
    # 8. 删除Be, H
    stacked_structure_withoutH = remove_hydrogen_and_convert(stacked_structure, element_order)
    write_lammps_dump_text("stacked_structure.lammpstrj", stacked_structure_withoutH)