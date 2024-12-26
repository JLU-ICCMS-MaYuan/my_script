#!/usr/bin/env python3

import argparse
from typing import Union, Sequence
import numpy as np

from ase.io import read, write
from ase.io.lammpsdata import write_lammps_data
from ase.calculators.lammps import Prism, convert
from ase import Atoms


def write_lammps_dump_text(
        filename: str,
        images: Union[Atoms, Sequence[Atoms]],
        specorder: list = None,
        write_image_flags: bool = False,
        reduce_cell: bool = False,
        units: str = 'metal',
        prismobj: Prism = None,
        ):
    symbols = images[0].get_chemical_symbols()
    totalatomsnum = len(images[0])
        
    with open(filename, 'w') as fd:
        for idx, image in enumerate(images):
            if specorder is None:
                # This way it is assured that LAMMPS atom types are always
                # assigned predictably according to the alphabetic order
                species = sorted(set(symbols))
            else:
                # To index elements in the LAMMPS data file
                # (indices must correspond to order in the potential file)
                species = specorder
                species2index = {element: idx for idx, element in enumerate(specorder)}
                sorted_indices = sorted(range(totalatomsnum), key=lambda i: species2index[image.get_chemical_symbols()[i]])
                image = image[sorted_indices]

            fd.write("ITEM: TIMESTEP\n")
            fd.write(f"{idx}\n")
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


def remove_hydrogen_and_convert(input_file, output_file, element_order):
    """
    读取 LAMMPS 轨迹文件，删除氢原子，按照指定顺序重新排列元素，写入 LAMMPS 格式文件。
    
    Args:
        input_file (str): 输入的 LAMMPS trajectory 文件名。
        output_file (str): 输出的 LAMMPS 格式文件名。
        element_order (list): 指定的元素排列顺序，例如 ['C', 'O', 'H', 'N']。
    """
    # 读取 LAMMPS 轨迹
    frames = read(input_file, format='lammps-dump-text', index=':')
    print(f"读取到 {len(frames)} 帧轨迹")

    # 处理每一帧，删除氢原子并按照指定顺序排列元素
    processed_frames = []
    for i, frame in enumerate(frames):
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

        processed_frames.append(new_frame)
        # print(f"帧 {i+1}: 保留 {len(new_frame)} 个原子")

    # write(output_file, processed_frames, format='vasp-xdatcar')
    
    # 将结果写入 LAMMPS 格式文件
    # write和write_lammps_data是等价的
    # write(output_file, processed_frames[0], format='lammps-data')
    # write_lammps_data(output_file, processed_frames, specorder=['La', 'Y', 'Ce', 'Th'])

    # write(output_file, processed_frames, format='lammps-dump-text') # 因为不存在 write_lammps_dump_text 这个函数，所以无法调用到 write_lammps_dump_text
    # 所以我自己写了一个
    write_lammps_dump_text(output_file, processed_frames)
    print(f"转换完成，输出文件为 {output_file}")


if __name__ == "__main__":


    # 设置输入输出文件
    parser = argparse.ArgumentParser(description="Process LAMMPS trajectory files.")

    # 添加输入文件和输出文件的参数
    parser.add_argument('-i', '--input_file', type=str, help="Input LAMMPS trajectory file")
    parser.add_argument('-o', '--output_file', type=str, help="Output LAMMPS trajectory file")

    # 解析命令行参数
    args = parser.parse_args()

    # 获取输入和输出文件路径
    input_file = args.input_file
    output_file = args.output_file

    # 定义元素排列顺序
    element_order = ['La', 'Y', 'Ce', 'Th']  # 按照这个顺序排列元素

    # 执行转换
    remove_hydrogen_and_convert(input_file, output_file, element_order)
