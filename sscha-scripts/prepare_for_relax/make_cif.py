from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter

# 读取 POSCAR 文件
poscar = Poscar.from_file("POSCAR")

# 获取结构对象
structure = poscar.structure

# 分析结构的空间群
spacegroup_analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = spacegroup_analyzer.get_symmetrized_structure()

# 将修正后的结构和对称性信息写入 CIF 文件
cif_writer = CifWriter(symmetrized_structure, symprec=0.01,write_magmoms=False)
cif_writer.write_file("structure.cif")

space_group_symbol = spacegroup_analyzer.get_space_group_symbol()

