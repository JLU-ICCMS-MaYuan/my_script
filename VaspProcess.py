#!/usr/bin/env python3
'''
Use examples:
    1. 分析结构优化的enthalpy
        VaspOutcar.py -d scripts_tests/201/50.0-75.0/ -opt
            -d scripts_tests/201/50.0-75.0/     -d 指定的目录下存放着大量的子目录, 
                                                这些子目录里面存放着每一个结构优化后的输出
            -opt                                   提取结构优化后的enthalpy
    2. 分析encut
        VaspOutcar.py -d encut/ -et -posf -plot
    3. 分析kspacing
        VaspOutcar.py -d kspacing/ -kp -posf -plot
'''
import os
import re
from collections import defaultdict
from argparse import ArgumentParser

from pymatgen.core.structure import Structure

def get_enthalpy(outcar_path):
    enthalpy_list = []
    with open(outcar_path, "r") as outcar:
        for line in outcar.readlines():
            if "enthalpy" in line:
                enthalpy = re.search(r"\-*\d+\.\d+", line).group()
                enthalpy_list.append(enthalpy)
    return enthalpy_list[-1]

def get_free_energy(outcar_path):
    free_energy_list = []
    with open(outcar_path, "r") as outcar:
        for line in outcar.readlines():
            if "free  energy   TOTEN" in line:
                free_energy = re.search(r"\-*\d+\.\d+", line).group()
                free_energy_list.append(free_energy)
    return free_energy_list[-1]    

def get_energy_without_entropy(outcar_path):
    energy_without_entropy_list = []
    with open(outcar_path, "r") as outcar:
        for line in outcar.readlines():
            if "energy  without entropy" in line:
                energy = re.search(r"\-*\d+\.\d+", line).group()
                energy_without_entropy_list.append(energy)
    return energy_without_entropy_list[-1]

def get_energy_sigma_0(outcar_path):
    energy_sigma_0_list = []
    with open(outcar_path, "r") as outcar:
        for line in outcar.readlines():
            if "energy(sigma->0)" in line:
                energy = re.search(r"\-*\d+\.\d+", line).group()
                energy_sigma_0_list.append(energy)
    return energy_sigma_0_list[-1]   

def plot(path_dat, path_jpeg):
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_table(path_dat, sep='\s+', header=0) # \s表示由空格作为分隔符, +表示有多个空格   index_col=0 是否指定第一列为索引
    # 若header = None，则表明数据中没有列名行；
    # 若header = 0，则表明第一行为列名
    df = df.sort_values(by=["name"])
    df['energy_per_atoms_sep']=df['energy_per_atoms'].diff()
    print(path_dat, "\n", df)
    df.to_csv(path_dat, index=False)
    # df.plot(x="name", y=["energy_per_atoms", "energy_per_atoms_sep"])
    df.plot(x="name", y=["energy_per_atoms"])
    plt.savefig(path_jpeg)

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--directory-vasp-run",
        default=None,
        type=str,
        # required=True,
        dest="directory_vasp_run",
        help="请输入结构所在的目录"
    )
    parser.add_argument(
        "-posf",
        "--poscar-file-read-flag",
        default=True,
        action="store_true",
        # required=True,
        dest="poscar_file_read_flag",
        help="是否读入目录中的POSCAR"
    )
    parser.add_argument(
        "-opt",
        "--structure-opt",
        default=False,
        action="store_true",
        dest="structure_opt",
        help="检查结构优化后的目录中的enthalpy"
    )
    parser.add_argument(
        "-encut",
        "--encut-test",
        default=False,
        action="store_true",
        dest="test_encut",
        help="检查测试encut后的能量"
    )
    parser.add_argument(
        "-kspacing",
        "--kspacing-test",
        default=False,
        action="store_true",
        dest="test_kspacing",
        help="检查测试kspacing后的能量"
    )
    parser.add_argument(
        "-sigma",
        "--sigma-test",
        default=False,
        action="store_true",
        dest="test_sigma",
        help="检查测试sigma后的能量"
    )
    parser.add_argument(
        '-plot',
        '--plot-picture-for-encuttest-kspacing-test',
        action='store_true',
        default=False,
        dest='plot_flag',
        help="是否为encut测试或者kspacing测试画图画图"
            )
    args = parser.parse_args()
    directory_vasp_run = args.directory_vasp_run
    structure_opt      = args.structure_opt
    test_encut         = args.test_encut
    test_kspacing      = args.test_kspacing
    test_sigma         = args.test_sigma
    plot_flag          = args.plot_flag

    if directory_vasp_run is not None:
        dir_vr = os.path.abspath(directory_vasp_run)
    
    if structure_opt:
        enthalpy_dict = defaultdict(dict)
        if os.path.exists("enthalpy.dat"):
            os.remove("enthalpy.dat")
        with open("enthalpy.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20} {:<20} \n".format("Number", "formula", "enthalpy", "enthalpy_per_atoms")) 

        i = 100000
        print("提取文件的路径如下：")
        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                print(root)
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                # print(outcar_path)
                struct = Structure.from_file(poscar_path)   
                atoms_amount = struct.composition.num_atoms
                formula = str(struct.composition.iupac_formula).replace(" ", "")
                comp_dict = dict(struct.composition.get_el_amt_dict())
                try:
                    enthalpy = get_enthalpy(outcar_path)
                    enthalpy_per_atoms = float(enthalpy) / atoms_amount
                    number = re.findall(r"\d{1,4}", root.split("/")[-1])
                    if number:
                        enthalpy_dict[formula]["Number"] = number[0]
                        number = number[0]
                    else:
                        i += 1
                        enthalpy_dict[formula]["Number"] = str(i) +"-"+ os.path.basename(root)
                        number =  str(i) +"-"+ os.path.basename(root)
                    enthalpy_dict[formula]["formula"]    = formula
                    enthalpy_dict[formula]["enthalpy"]   = enthalpy
                    enthalpy_dict[formula]["enthalpy/atoms"] = enthalpy_per_atoms
                    for key, value in comp_dict.items():
                        enthalpy_dict[formula][key] = value
                except Exception as e:
                    print(e.args)
                    print(f"error in {root}")
                    continue                

        import pandas as pd
        # import pprint
        # pprint.pprint(enthalpy_dict); input()
        df = pd.DataFrame(enthalpy_dict).T  # 分隔符的用法: \s表示由空格作为分隔符, +表示有多个空格
        # 将所有的 NaN 值替换为 0
        df = df.fillna(0)

        try:
            print("尝试处理端点值")
            # (df.iloc[:, 4:] != 0) 可以返回一个布尔值的df，非零的位置是True，零的位置是False
            # sum(axis=1) 对每一行相加（例如某一行是 False，True，True对应着0, 18, 0）,那么相加的结构是1就说明只有一个元素，说明提取到的这一行是单质的行。
            endpoints_df = df[(df.iloc[:, 4:] != 0).sum(axis=1) == 1]

            # 保留端点值中能量最低的结构，其余结构和能量都删除
            # df['formula'].str.extract('([A-Za-z]+)', expand=False) 该行代码可以用于提取化学式中的化学元素符号。 expand=False就是返回第一个匹配到的内容， expand=True是返回所有匹配到的内容
            idx = endpoints_df.groupby(endpoints_df['formula'].str.extract('([A-Za-z]+)', expand=False))['enthalpy/atoms'].idxmin()
            endpoints_df = endpoints_df.loc[idx]
            # 处理化合物的点
            comppoints_df = df[(df.iloc[:, 4:] != 0).sum(axis=1) > 1]

            print("端点处理好了，显示如下：")
            print(endpoints_df)

            # 拼接两个pandas
            df = pd.concat([endpoints_df, comppoints_df], axis=0)
        except Exception as e:
            print("端点值处理失败，可能你提供的结构内不包含端点值")
            print(e.args)
  
        df.sort_values(
            by=["formula", "enthalpy/atoms"],
            ascending=[True, True], # 按照升序排列
            inplace=True
            )
        # index=False表示不输出行索引 
        # header=False表示不输出列名
        # line_terminator='\n'表示使用换行符作为行分隔符
        # 最后一个参数justify='left'指定左对齐。
        df_string = df.to_string(index=False, justify='right', col_space=8)
        # 将字符串写入文件
        with open('enthalpy.dat', 'w') as file:
            file.write(df_string)
        # df.to_csv("enthalpy_sorted.csv", sep=' ', index=False, header=False, line_terminator='\n', justify='left')

    if test_encut:
        if os.path.exists("encut.dat"):
            os.remove("encut.dat")        
        with open("encut.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20}\n".format("name", "energy", "energy_per_atoms"))     
        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")

                struct = Structure.from_file(poscar_path);  atoms_amount     = struct.composition.num_atoms
                try:
                    energy = get_free_energy(outcar_path)
                    energy_per_atoms = float(energy) / atoms_amount
                except:
                    print(root)
                    continue
                name   = (root.split("/")[-1])
                with open("encut.dat", "a") as file:
                    file.write("{:<10} {:<20} {:<20}\n".format(name, energy, energy_per_atoms))
    if test_kspacing:
        if os.path.exists("kspacing.dat"):
            os.remove("kspacing.dat")
        with open("kspacing.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20}\n".format("name", "energy", "energy_per_atoms"))    
        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                print(outcar_path)
                struct = Structure.from_file(poscar_path);  atoms_amount     = struct.composition.num_atoms
                energy = get_free_energy(outcar_path)    ;  energy_per_atoms = float(energy) / atoms_amount
                name   = root.split("/")[-1]
                with open("kspacing.dat", "a") as file:
                    file.write("{:<10} {:<20} {:<20}\n".format(name, energy, energy_per_atoms))      
    if test_sigma:
        if os.path.exists("sigma.dat"):
            os.remove("sigma.dat")
        with open("sigma.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20}\n".format("name", "energy", "energy_per_atoms"))    
        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                print(outcar_path)
                struct = Structure.from_file(poscar_path);  atoms_amount     = struct.composition.num_atoms
                energy = get_free_energy(outcar_path)    ;  energy_per_atoms = float(energy) / atoms_amount
                name   = root.split("/")[-1]
                with open("sigma.dat", "a") as file:
                    file.write("{:<10} {:<20} {:<20}\n".format(name, energy, energy_per_atoms))      

    if os.path.exists("encut.dat") and plot_flag and test_encut:
        plot("encut.dat", "encut.jpeg")
    elif os.path.exists("kspacing.dat") and plot_flag and test_kspacing:
        plot("kspacing.dat", "kspacing.jpeg")
    elif os.path.exists("sigma.dat") and plot_flag and test_sigma:
        plot("sigma.dat", "sigma.jpeg")