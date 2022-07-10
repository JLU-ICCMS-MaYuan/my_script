#!/work/home/mayuan/miniconda3/envs/pyxtal/bin/python3

from pymatgen.core.structure import Molecule
from pymatgen.symmetry.groups import PointGroup
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

if __name__ == "__main__":
    pg = PointGroup("-43m")
    print(pg)
    
