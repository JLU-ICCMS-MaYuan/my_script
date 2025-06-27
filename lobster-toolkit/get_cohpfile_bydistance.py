#!/usr/bin/env python
import os
import sys
import argparse
from pymatgen.core.structure import Structure

# 设置命令行参数解析器
def parse_args():
    parser = argparse.ArgumentParser(description="生成LOBSTER输入文件")
    parser.add_argument('-e', '--COHPEnergy', nargs="+", type=float, help="COHP起始终止能量")
    parser.add_argument('-s', '--species_custom', nargs="+", type=str, help="元素1, 2")
    parser.add_argument('-z', '--zvalances', nargs="+", type=str, help="元素1, 2 及其价电子轨道, 例如 -z Nb:4s,4p,4d,5s H:1s")
    parser.add_argument('-d', '--d_limit', nargs="+", type=float, help="原子间最小最大距离")
    parser.add_argument('-m', '--mode', type=int, required=True, help="选择生成的输入文件版本: 5 或 0")
    return parser.parse_args()

# 解析价电子轨道参数 '-z'
def parse_zvalances(zvalances):
    zval_dict = {}
    for entry in zvalances:
        element, orbitals = entry.split(":")
        zval_dict[element] = ' '.join(orbitals.split(","))
    return zval_dict

# 写入lobsterin文件
def write_lobsterin(dirname, mode=5, COHPstartEnergy=None, COHPendEnergy=None, species_custom1=None, species_custom2=None, lower_d=None, upper_d=None, struct=None, zval_dict=None):
    if mode == 5:
        lobsterin_path = os.path.join(dirname, "lobsterin")
        with open(lobsterin_path, "w") as f:
            f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy))
            f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
            f.write('usebasisset pbeVaspFit2015\n')  # 基组
            for spe in struct.types_of_specie:
                f.write('basisfunctions {} {} \n'.format(spe.name, zval_dict[spe.name]))   # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
            # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
            # lobster认为你指定的原子对是不具有周期性的
            # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
            # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
            # 所以这里我抛弃了设置原子对来计算成键强度的方法。
            # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
            # for pair, idx, d in pairs_idxs_d:
            #     if lower_d <= d <= upper_d:
            #         f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
            f.write("cohpGenerator from {} to {} type {} type {} orbitalWise\n".format(lower_d, upper_d, species_custom1, species_custom2))
    elif mode == 0:
        lobsterin_path = os.path.join(dirname, "lobsterin")
        with open(lobsterin_path, "w") as f:
            f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
            f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
            f.write('usebasisset pbeVaspFit2015\n')  #  # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
            f.write('gaussianSmearingWidth 0.05\n')
            for spe in struct.types_of_specie:
                f.write('basisfunctions {} {}\n'.format(spe.name, zval_dict[spe.name]))   # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
            # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
            # lobster认为你指定的原子对是不具有周期性的
            # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
            # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
            # 所以这里我抛弃了设置原子对来计算成键强度的方法。
            # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
            # for pair, idx, d in pairs_idxs_d:
            #     if lower_d <= d <= upper_d:
            #         f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
            f.write("cohpGenerator from {} to {} type {} type {} orbitalWise\n".format(lower_d, upper_d, species_custom1, species_custom2))
    else:
        print(f"mode = {mode}, The scripts doesn't support it")
# 主函数
def main():
    # 获取命令行参数
    args = parse_args()

    # 提取参数
    COHPstartEnergy = args.COHPEnergy[0]
    COHPendEnergy = args.COHPEnergy[1]
    species_custom1 = args.species_custom[0]
    species_custom2 = args.species_custom[1]
    lower_d = args.d_limit[0]
    upper_d = args.d_limit[1]
    mode = args.mode  # 选择生成的输入文件版本

    # 解析价电子轨道
    zval_custom = parse_zvalances(args.zvalances)

    # 输出解析后的价电子信息（可选）
    print("Parsed Z-Valences:")
    for element, orbitals in zval_custom.items():
        print(f"{element}: {orbitals}")

    # 读取结构文件
    struct = Structure.from_file("POSCAR")

    # 创建输出目录
    dirs = "{}_{}_{}_{}".format(species_custom1, species_custom2, lower_d, upper_d)
    if not os.path.exists(dirs):
        os.mkdir(dirs)

    # 创建符号链接
    files = ['WAVECAR', 'CONTCAR', 'KPOINTS', 'OUTCAR', 'POTCAR', 'vasprun.xml']
    for file in files:
        if os.path.exists(file):
            os.system(f"ln -s {os.path.abspath(file)} {dirs}")
        else:
            print(f"{file} does not exist.")

    # 写入LOBSTER输入文件
    write_lobsterin(dirs, mode=mode, COHPstartEnergy=COHPstartEnergy, COHPendEnergy=COHPendEnergy, species_custom1=species_custom1, species_custom2=species_custom2, lower_d=lower_d, upper_d=upper_d, struct=struct, zval_dict=zval_custom)

    # 输出信息
    info = """Note: -------------------------------
1. 如果运行了lobster得不到COHP，说明NBANDS不够多
2. 如果运行了vasp发现不收敛，可以尝试做如下改变
    ISTART  = 0   ---->  ISTART  = 1   读取波函数
    ICHARG  = 2   ---->  ICHARG  = 11  读取电荷密度
"""
    print(info)

# 确保仅在脚本作为主程序运行时才执行
if __name__ == "__main__":
    main()

