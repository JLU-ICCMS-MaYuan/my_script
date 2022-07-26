#!/work/home/mayuan/miniconda3/envs/pyxtal/bin/python3
'''
The requirement is:
    1. 赝势文件UPFs收集好之后放在一个子目录中, 例如放在名为pp的子目录中
    2. 当前目录要放置好POSCAR结构。目前的脚本还只能读取POSCAR,
    3. 读取POSCAR后程序会自动将POSCAR变为原胞, 然后输出原胞的PPOSCAR和晶胞的BPOSCAR
    4. I have no idea at the moment
The Program Logic:
    1. 在当前目录下准备好POSCAR
    2. 读取POSCAR, 读取赝势目录
    3. 写qe进行单点能计算的输入文件scf.in
Use examples:
    vasp2qe_relax.py -ppd pp -p 1000
        -ppd pp         指定赝势目录
'''

import os
import re
from argparse import ArgumentParser

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-ppd",
        "--pseudopotential-directory",
        type=str,
        dest="pp_dir",
        required=True,
        help="给出赝势存放的路径"
    )
    args = parser.parse_args()
    pp_dir = args.pp_dir
    
    pp_dir_path = os.path.abspath(pp_dir)
    for root, dirs, files in os.walk(pp_dir_path):
        pp_file = files


    # 读入结构文件POSCAR，导出PPOSCAR，BPOSCAR
    poscar_file = os.path.abspath("POSCAR")
    if not os.path.exists(poscar_file):
        raise FileExistsError("POSCAR doesn't exist!")
    struct = Structure.from_file("POSCAR")
    spa = SpacegroupAnalyzer(struct)
    bstruct = spa.get_conventional_standard_structure()
    pstruct = spa.get_primitive_standard_structure()
    Poscar(bstruct).write_file("BPOSCAR")
    Poscar(pstruct).write_file("PPOSCAR")


    # 处理PPOSCAR的pymatgen对象
    # 获得元素名称 和 每种元素的原子个数
    composition = pstruct.composition.get_el_amt_dict()
    # 获得体系的 化学配比
    system_name = pstruct.composition.formula
    # 获得元素种类的个数
    species_quantity = len(composition)
    # 获得每种元素的相对原子质量
    all_atoms_quantity = int(sum(composition.values()))
    # 获得晶格矩阵
    cell_parameters    = pstruct.lattice.matrix
    # 获得原子分数坐标
    fractional_sites  = pstruct.sites
    with open("scf.in", "w") as qe:
        qe.write("&CONTROL\n")
        qe.write(" calculation='scf',              \n")
        qe.write(" restart_mode='from_scratch',    \n")
        qe.write(" prefix='{}',                    \n".format(system_name.replace(' ', '')))
        qe.write(" pseudo_dir='{}',                \n".format(pp_dir_path))
        qe.write(" outdir='./tmp',                 \n")
        qe.write(" forc_conv_thr = 1.0d-3,         \n")
        qe.write(" etot_conv_thr = 1.0d-4,         \n")
        qe.write("/\n")

        qe.write("&SYSTEM\n")
        qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
        qe.write(" nat={},                         \n".format(all_atoms_quantity))
        qe.write(" ntyp={},                        \n".format(species_quantity))
        qe.write(" occupations = 'smearing',       \n")
        qe.write(" smearing = 'methfessel-paxton'  \n")
        qe.write(" degauss = 0.02                  \n")
        qe.write(" ecutwfc = 60,                   \n")
        qe.write(" ecutrho = 720,                  \n")
        qe.write("/\n")

        qe.write("&ELECTRONS\n")
        qe.write(" conv_thr = 1.0d-9               \n")
        qe.write(" mixing_beta = 0.8d0             \n")
        qe.write("/\n")

        qe.write("ATOMIC_SPECIES                   \n")
        for species_name in composition.keys():
            for species_pseudo in pp_file:
                if species_name.lower() in species_pseudo.lower():
                    element      = Element(species_name)
                    species_mass = str(element.atomic_mass).strip("amu")
                    qe.write(" {:<5}  {:<10}  {:<50} \n".format(species_name, species_mass, species_pseudo))
        qe.write("CELL_PARAMETERS(angstrom)        \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
        for cell_p in cell_parameters:
            cell_p = list(map(str, cell_p))
            qe.write("{:>25} {:>25} {:>25} \n".format(cell_p[0], cell_p[1], cell_p[2]))
        qe.write("ATOMIC_POSITIONS (crystal)       \n")
        for site in fractional_sites:
            coord = list(map(str, site.frac_coords))
            name  = re.search(r"[A-Za-z]+", str(site.species)).group()
            # 左对齐5个字符，左对齐30个字符
            qe.write("{:<5} {:<30} {:<30} {:<30} \n".format(name, coord[0], coord[1], coord[2]))

