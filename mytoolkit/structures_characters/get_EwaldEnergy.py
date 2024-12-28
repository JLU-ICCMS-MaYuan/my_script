#!/usr/bin/env python3
import os
import argparse
from pymatgen.core import Structure
from pymatgen.analysis.energy_models import EwaldSummation

# 设置命令行参数解析器
def parse_arguments():
    parser = argparse.ArgumentParser(description="Get file paths based on input directory and conditions.")
    parser.add_argument('-i', '--input-file-or-directory', type=str, help="The structure file (e.g., POSCAR, CONTCAR) or a directory")
    parser.add_argument('-c', '--detailed-conditions', type=str, nargs='+', default=[], 
                        help="List of subdirectories or patterns for detailed selection")
    parser.add_argument('-d', '--oxidation_data', type=str, nargs='+', help="List of elements and their oxidation states (e.g., -d Nb:4 H:-1)")
    return parser.parse_args()

def parse_data(data_list):
    # 将输入的data转换为字典形式
    oxidation_data = {}
    for item in data_list:
        try:
            element, oxidation_state = item.split(':')
            oxidation_data[element] = int(oxidation_state)  # 将氧化态转为整数
        except ValueError:
            print(f"Invalid format for element and oxidation state: {item}. It should be 'Element:oxidation_state'")
    return oxidation_data

# 设置需要计算的文件类型，假设所有结构文件都为 .vasp 文件
def calculate_electric_energy(structure_file, oxidation_data):
    try:
        # 从文件加载结构
        structure = Structure.from_file(structure_file)
        structure.add_oxidation_state_by_element(oxidation_data)
        # 计算静电能
        energy = EwaldSummation(structure).total_energy
        
        # 返回计算结果
        return energy
    except Exception as e:
        print(f"Error in processing {structure_file}: {e}")
        return None

def get_files(input_file_or_directory, detailed_conditions):
    file_paths = []

    # 检查输入的是文件还是目录
    if os.path.isdir(input_file_or_directory):
        # 如果是目录，遍历所有文件
        for root, dirs, files in os.walk(input_file_or_directory):
            for filename in files:
                filepath = os.path.join(root, filename)
                # 如果文件路径符合所有条件
                if all(pattern in filepath for pattern in detailed_conditions):
                    file_paths.append(filepath)
    elif os.path.isfile(input_file_or_directory):
        # 如果是单个文件
        file_paths.append(input_file_or_directory)
    else:
        print(f"The path {input_file_or_directory} is neither a valid file nor a directory.")

    return file_paths

def main():
    # 解析命令行参数
    args = parse_arguments()
    input_file_or_directory = args.input_file_or_directory
    detailed_conditions = args.detailed_conditions  # 获取详细的筛选条件

    # 解析 oxidation_data 字典
    oxidation_data = parse_data(args.oxidation_data)

    # 获取符合条件的文件路径
    files = get_files(input_file_or_directory, detailed_conditions)

    # 用来存储文件和计算出来的静电能
    energy_results = []

    # 遍历所有文件并计算静电能
    for file in files:
        energy = calculate_electric_energy(file, oxidation_data)
        if energy is not None:
            energy_results.append((file, energy))

    # 按照静电能从小到大排序
    sorted_energy_results = sorted(energy_results, key=lambda x: x[1])

    # 输出结果
    print(f"file Electrostatic Energy(eV)")
    for file, energy in sorted_energy_results:
        print(f"{file:<}   {energy:>.6f} eV")

if __name__ == '__main__':
    main()
