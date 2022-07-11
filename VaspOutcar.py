#!/work/home/may/miniconda3/bin/python3
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

def plot_encut_or_kspacing(path_dat, path_jpeg):
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_table(path_dat, sep='\s+', header=0) # \s表示由空格作为分隔符, +表示有多个空格   index_col=0 是否指定第一列为索引
    # 若header = None，则表明数据中没有列名行；
    # 若header = 0，则表明第一行为列名
    df = df.sort_values(by=["name"])
    df['energy_per_atoms_sep']=df['energy_per_atoms'].diff()
    print(path_dat, "\n", df)
    df.plot(x="name", y=["energy_per_atoms", "energy_per_atoms_sep"])
    
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
        "-et",
        "--encut-test",
        default=False,
        action="store_true",
        dest="test_encut",
        help="检查测试encut后的能量"
    )
    parser.add_argument(
        "-kp",
        "--kspacing-test",
        default=False,
        action="store_true",
        dest="test_kspacing",
        help="检查测试kspacing后的能量"
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
    plot_flag          = args.plot_flag

    if directory_vasp_run is not None:
        dir_vr = os.path.abspath(directory_vasp_run)
    
    if structure_opt:
        enthalpy_dict = defaultdict(dict)
        if os.path.exists("enthalpy.dat"):
            os.remove("enthalpy.dat")
        with open("enthalpy.dat", "w") as file:
            file.write("{:<60} {:<20} {:<20} \n".format("formula", "enthalpy", "enthalpy_per_atoms")) 

        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                # print(outcar_path)
                struct   = Structure.from_file(poscar_path);    atoms_amount       = struct.composition.num_atoms;
                enthalpy = get_enthalpy(outcar_path)       ;    enthalpy_per_atoms = float(enthalpy)/atoms_amount
                name     =  root.split("/")[-1]            ;    formula  = re.search(r"[A-Za-z]{1,2}\d+[A-Za-z]{1,2}\d+", name).group()

                enthalpy_dict[name]["formula"]            = formula
                enthalpy_dict[name]["enthalpy"]           = enthalpy
                enthalpy_dict[name]["enthalpy_per_atoms"] = enthalpy_per_atoms
                # with open("enthalpy.dat", "a") as file:
                #     file.write("{:<60} {:<20} {:<20} \n".format(root.split("/")[-1], enthalpy, enthalpy_per_atoms)) 
                
        import pandas as pd
        # import pprint
        # pprint.pprint(enthalpy_dict); input()
        df = pd.DataFrame(enthalpy_dict).T  # 分隔符的用法: \s表示由空格作为分隔符, +表示有多个空格
        df.sort_values(
            by=["enthalpy"],
            ascending=[True], # 按照升序排列
            inplace=True
            )
        df.to_excel("enthalpy_sorted.xlsx")

    if test_encut:
        if os.path.exists("test_encut.dat"):
            os.remove("test_encut.dat")
        with open("test_encut.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20}\n".format("name", "energy", "energy_per_atoms"))
        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                print(outcar_path)
                struct = Structure.from_file(poscar_path);  atoms_amount     = struct.composition.num_atoms
                energy = get_free_energy(outcar_path)    ;  energy_per_atoms = float(energy) / atoms_amount
                name   = (root.split("/")[-1])
                with open("test_encut.dat", "a") as file:
                    file.write("{:<10} {:<20} {:<20}\n".format(name, energy, energy_per_atoms))

    if test_kspacing:
        if os.path.exists("test_kspacing.dat"):
            os.remove("test_kspacing.dat")
        with open("test_kspacing.dat", "w") as file:
            file.write("{:<10} {:<20} {:<20}\n".format("name", "energy", "energy_per_atoms"))    

        for root, dirs, files in os.walk(dir_vr):
            if "OUTCAR" in files and "POSCAR" in files:
                outcar_path = os.path.join(root, "OUTCAR")
                poscar_path = os.path.join(root, "POSCAR")
                print(outcar_path)
                struct = Structure.from_file(poscar_path);  atoms_amount     = struct.composition.num_atoms
                energy = get_free_energy(outcar_path)    ;  energy_per_atoms = float(energy) / atoms_amount
                name   = root.split("/")[-1]
                with open("test_kspacing.dat", "a") as file:
                    file.write("{:<10} {:<20} {:<20}\n".format(name, energy, energy_per_atoms))      


    if os.path.exists("test_encut.dat") and plot_flag and test_encut:
        plot_encut_or_kspacing("test_encut.dat", "test_encut.jpeg")
    elif os.path.exists("test_kspacing.dat") and plot_flag and test_kspacing:
        plot_encut_or_kspacing("test_kspacing.dat", "test_kspacing.jpeg")