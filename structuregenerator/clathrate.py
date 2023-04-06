#!/usr/bin/env python
# python >= 3.9

import sys
from collections import defaultdict
from itertools import chain, groupby
from typing import NoReturn, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.analysis import local_env
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class HydrideClathrate:
    STRATEGY = local_env.CrystalNN(
        distance_cutoffs=None,
        x_diff_weight=-1,
        porous_adjustment=False,
    )
    def __init__(
        self,
        structure: Structure,
        min_cn: int = 4,
    ):
        """A class to determine a structure is a cathrate or not

        1. split atoms by hydrogen - [H], and non-hydrogen - [core]
        2. find the cages around [core] by Voronoi algorithm (CrystalNN)
        3. whether there remains a H which do not form a cage
        4. whether each cage vertex is shared by different cages
            - shr_num == 2  : most likely a covalent network structure
            - shr_num >= 3  : most likely a staggered stacked

        key attributes:
            - remain_H, remain_H_ratio  -- index of which H do not form cage
            - hybrid  -- Fale if all cage vertexes are H, else True
            - shr_num, shr_num_avg  -- share number of each vertex
            - describe()

        Note: we ignore the specie of cage vertex when finding cages.

        Args:
            structure(Sturcutre): pymatgen structure.Structure
            min_cn(int): minimal coordination number, i.e. number of cage
            vertex, default 4
        """
        self.structure = structure
        self.min_cn = min_cn

        self.idx_split = self._split_index_by_H(self.structure)
        self._check_hydride()  # check hydride
        self._check_dist(0.85, 1.9)
        self.ions_raddi = {
            "Ca":1.14,
            "Mg":0.86,
            "La":1.17,
            "Ce":1.15,
            "Li":0.9,
            "Y":1.04,
            "Th":1.08,
            # "Be":
            # "B":
        }
    @classmethod
    def from_file(cls, file_path):
        structure = Structure.from_file(file_path)
        self = cls(structure)
        return self

    def __repr__(self):
        return self.describe().__repr__()

    @staticmethod
    def _split_index_by_H(structure) -> dict[str, list[int]]:
        """split index by is H or not

        ```text
        {'H': [...], 'core': [...]}
               ^indexes       ^indexes
        ```
        """
        def keyfunc(k):
            if k[1].specie.symbol == 'H':
                return 'H'
            else:
                return 'core'

        # 将keyfunc函数作用于序列的各个元素(这里的各个元素就是指Structure结构里面各个原子的坐标构成的类PeriodicSites class)
        # 根据keyfunc函数结果(它只有两种结果'H'和'core')
        # 将拥有相同函数结果的元素分到一个新的迭代器(将返回结果为'H'的PeriodicSites类放入同一个迭代器，结果为'core'的PeriodicSites类的放入另一个迭代器)。
        # 每个新的迭代器以函数返回结果为标签。(针对当前keyfunc函数，只有两种类型的迭代器，其标签分别为'H'和'core')
        return {
            k: [i[0] for i in g]
            for k, g in groupby(enumerate(structure), keyfunc) # enumerate(structure)的输出结果大概是这样的：(序号-从0开始, PeriodicSite: H (2.8002, 0.0000, 0.0000) [0.7260, 0.0000, 0.0000])
        }

    def _check_hydride(self):
        """check if the structure is a hydride or not"""
        if len(self.idx_split) != 2:
            self.hydride = False
            print(f"{self.formula} is Not a hydride")
        else:
            self.hydride = True

    def _check_dist(self, h2h_dist_lowlimit=0.85, nonh2nonh_dist_lowlimit=1.9):
        '''
        1. check h2h_1nearest > 0.8 
        2. check nonh2nonh_1nearest > 2.0
        '''
        if (self.h2h_1nearest[0][0] <= h2h_dist_lowlimit) or (self.nonh2nonh_1nearest[0][0] <= nonh2nonh_dist_lowlimit):
            self.check_dist = False
        else:
            self.check_dist = True
        
        
    # basic properties
    @property
    def spg_symmstry(self):
        spginfo    = SpacegroupAnalyzer(self.structure)
        spg_num    = spginfo.get_space_group_number()
        spg_symbol = spginfo.get_space_group_symbol()
        self._spg_symmstry = spg_symbol+"("+str(spg_num)+")"
        return self._spg_symmstry

    @property
    def h_structure(self):
        self._h_structure = self.structure.copy()
        self._h_structure.remove_sites(self.idx_split['core'])
        return self._h_structure

    @property
    def nonh_structure(self):
        self._nonh_structure = self.structure.copy()
        self._nonh_structure.remove_sites(self.idx_split['H'])
        return self._nonh_structure

    @property
    def formula(self):
        self._formule = self.structure.composition.reduced_formula # 本质上reduced_composition调用的还是reduced_formula
        return self._formule

    @property
    def graph(self) -> np.ndarray:
        """graph in np.ndarray

        [i, j, jix, jiy, jiz]
        (Nnode * 5)

        Note: g.edges(data='to_jimage') returns [(i, j, (to_jimage)), ...]

        Note: graph of reversed direction is concatenated

        """
        if not hasattr(self, '_graph'):
            g = StructureGraph.with_local_env_strategy(
                self.structure, self.__class__.STRATEGY).graph
            edges = [[i, j] + list(to_jimage)
                     for i, j, to_jimage in g.edges(data='to_jimage')]
            edges = np.array(edges, dtype=int)
            # reverse edge direction and flip to_jimage
            rev_edges = (
                edges[:, [1, 0, 2, 3, 4]] * np.array([[1, 1, -1, -1, -1]])
            )
            self._graph = np.concatenate([edges, rev_edges])
        return self._graph

    @property
    def cages(self) -> dict:
        """cages of each core atom

        Return:
            {core_idx: <cage>}, where <cage>={'vertex': N, 'vector': N*3}
        """
        if not hasattr(self, '_cages'):
            cages = [
                self.get_cage(core_idx) for core_idx in self.idx_split['core']
            ]
            self._cages = {
                core_idx: cage
                for core_idx, cage in zip(self.idx_split['core'], cages)
                if cage is not None
            }
        return self._cages

    def get_cage(self, idx) -> Optional[dict]:
        """Give one atom index as core, get the cage index around it

        Return vertex index of the cage, and vector from core to each vertex.

        None: if the number of vertex atoms is less then <self.min_cn>,
        a cage cannot be formed or we do not treat it as a cage.

        Return:
            {'vertex': N, 'vector': N*3}: if form a cage
            None: if not form a cage
        """
        cage = self.graph[np.asarray(self.graph[:, 0] == idx).nonzero()]
        if len(cage) < self.min_cn:
            # print(f"Not a cage around index: {idx}", file=sys.stderr)
            return None
        i, j = cage[:, 0], cage[:, 1]
        to_jimage = cage[:, -3:]
        v = - self.structure.cart_coords[i] \
            + self.structure.cart_coords[j] \
            + np.einsum('ni,ij->nj', to_jimage, self.structure.lattice.matrix)
        return {'vertex': j, 'vector': v, 'to_jimage': to_jimage}

    def get_cage_mol(self, idx)-> Molecule:
        hcage_coord = self.cages[idx]['vector']
        species = [self.structure[self.cages[idx]['vertex'][_]].specie.name for _ in range(len(hcage_coord))] + [self.structure[idx].specie.name]
        coords  = np.concatenate((hcage_coord, np.array([[0., 0., 0.]])), axis=0)
        cage_mol= Molecule(species, coords)
        return cage_mol


    # Cages properties
    @property
    def remain_H(self) -> set:
        """H index of which not form a cage

        Returns:
            set: H index
        """
        if not hasattr(self, '_remain_H'):
            all_cages_vertex = chain.from_iterable(
                cage['vertex'] for cage in self.cages.values())
            self._remain_H = set(self.idx_split['H']) - set(all_cages_vertex)
        return self._remain_H

    @property
    def remain_H_ratio(self) -> float:
        return len(self.remain_H) / len(self.idx_split['H'])

    @property
    def hybrid_cage(self):
        """any cage whose vertexes are not all H

        Returns:
            {cage_idx: [non-hydrogen-vertex-idx]}
        """
        if not hasattr(self, '_hybrid_cage'):
            self._hybrid_cage = {}
            for core_idx, cage in self.cages.items():
                # non-hydrogen-vertex-idx
                v = set(cage['vertex']) - set(self.idx_split['H'])
                if len(v) > 0:
                    self._hybrid_cage[core_idx] = list(v)
        return self._hybrid_cage

    @property
    def hybrid(self) -> bool:
        if len(self.hybrid_cage) == 0:
            return False
        else:
            return True

    @property
    def shr_num(self) -> dict:
        """share number of each vertex"""
        if not hasattr(self, '_shr_num'):
            all_cages_vertex = chain.from_iterable(
                cage['vertex'] for cage in self.cages.values())
            self._shr_num = defaultdict(int)
            for vertex in all_cages_vertex:
                self._shr_num[vertex] += 1
        return self._shr_num

    @property
    def shr_num_avg(self) -> np.ndarray:
        """average share number"""
        return np.mean(list(self.shr_num.values()))

    @property
    def cage_regularity(self):
        """Check the regularity of the cages in a hydride"""
        if not hasattr(self, '_regularity'):
            self._cage_regularity = defaultdict(float)
            for core_idx, cage in self.cages.items():
                dist_core2h = np.linalg.norm(cage['vector'], axis=1)
                self._cage_regularity[core_idx] = np.std(dist_core2h)
        return self._cage_regularity

    @property
    def cage_regularity_avg(self):
        """average regularity"""
        return np.mean(list(self.cage_regularity.values()))


    # H subcrystal properties
    @property
    def h_network(self) -> Optional[dict]:
        '''The distance between H-H is counted'''
        if not hasattr(self, '_h_network'):
            _h_network = [
                self.get_neighbers(center_h_vertex) for center_h_vertex in self.idx_split['H']
            ]
            self._h_network = {
                center_h_vertex: neighbers
                for center_h_vertex, neighbers in zip(self.idx_split['H'], _h_network)
                if neighbers is not None
        }
        return self._h_network

    def get_neighbers(self, idx):
        '''Get all the hydrogen atoms around an h atom index `idx`'''
        neighbers = self.graph[ 
            np.asarray(
                (self.graph[:, 0] == idx) & (~np.isin(self.graph[:, 1], self.idx_split['core']))
                ).nonzero()
            ]
        i, j = neighbers[:, 0], neighbers[:, 1]
        to_jimage = neighbers[:, -3:]
        v = - self.structure.cart_coords[i] \
            + self.structure.cart_coords[j] \
            + np.einsum('ni,ij->nj', to_jimage, self.structure.lattice.matrix)
        return {'vertex': j, 'vector': v, 'to_jimage': to_jimage}

    @property
    def h2h_dists(self):
        """
        All H-H distances were obtained
        return self._h2h_dists, its shape is N (1D ndarray)
        
        """
        self._h2h_dists = self.h_structure.distance_matrix
        self._h2h_dists = self._h2h_dists[np.triu_indices(self._h2h_dists.shape[0], k=1)]
        self._h2h_dists =  np.sort(self._h2h_dists)
        return self._h2h_dists
    
    @property
    def h2h_1nearest(self):
        return self.Nnearest(1, self.h2h_dists, rtol=0.01)

    @property
    def h2h_2nearest(self):
        return self.Nnearest(2, self.h2h_dists, rtol=0.01)

    @property
    def h2h_3nearest(self):
        return self.Nnearest(3, self.h2h_dists, rtol=0.01)

    @property
    def h2h_4nearest(self):
        return self.Nnearest(4, self.h2h_dists, rtol=0.01)

    @property
    def h2h_5nearest(self):
        return self.Nnearest(5, self.h2h_dists, rtol=0.01)

    @property
    def h2h_network_regularity(self):
        if not hasattr(self, '_h2h_network_regularity'):
            self._h2h_network_regularity = defaultdict(float)
            for center_h_vertex, neighbers in self.h_network.items():
                dist_core2h = np.linalg.norm(neighbers['vector'], axis=1)
                self._h2h_network_regularity[center_h_vertex] = np.std(dist_core2h)
        return self._h2h_network_regularity

    @property
    def h2h_network_regularity_avg(self):
        """average regularity"""
        return np.mean(list(self.h2h_network_regularity.values()))
    
    
    # nonh subcrystal properties
    @property
    def nonh2nonh_dist(self):
        """
        All nonh2nonh distances were obtained
        return self._nonh2nonh_dist, its shape is N (1D ndarray)
        """
        self._nonh2nonh_dist = self.nonh_structure.distance_matrix
        self._nonh2nonh_dist = self._nonh2nonh_dist[np.triu_indices(self._nonh2nonh_dist.shape[0], k=1)]
        self._nonh2nonh_dist =  np.sort(self._nonh2nonh_dist)
        if self._nonh2nonh_dist.size == 0:
            self._nonh2nonh_dist = np.sort(np.array([self.nonh_structure.lattice.a, self.nonh_structure.lattice.b, self.nonh_structure.lattice.c]))
            return self._nonh2nonh_dist
        else:
            return self._nonh2nonh_dist

    @property
    def nonh2nonh_1nearest(self):
        return self.Nnearest(1, self.nonh2nonh_dist, rtol=0.01)

    @property
    def nonh2nonh_2nearest(self):
        return self.Nnearest(2, self.nonh2nonh_dist, rtol=0.01)

    @staticmethod
    def Nnearest(N, dist, rtol=0.1):
        '''
        Calculate H-H nearest distance of N order.
            if N=1, output the first  nearest H-H distance
            if N=2, output the second nearest H-H distance
            if N=3, output the third  nearest H-H distance

            if there is no existance of fifth order nearest distance, it will return  [[10000], 1]
            ...
        '''
        _dists = np.copy(dist)
        i = 0
        while True:
            diff = _dists - np.min(_dists)
            N_nearest = dist[(diff<rtol)&(diff>=0)]
            i += 1
            if   i == 1 and   np.all((diff<rtol)&(diff>=0)):
                return [N_nearest, len(N_nearest)]
            elif i >= N and (~np.all((diff<rtol)&(diff>=0))):
                return [N_nearest, len(N_nearest)]
            elif i >= N and i > 1 and np.all((diff<rtol)&(diff>=0)):
                return [[10000.0], 1]
            else:
                _dists[diff<rtol] = 10000.0
     
    @property
    def fraction_of_hydrogen_volume(self):
        """
        calculate hydrogen atoms 
        """
        ele_amt = self.structure.composition.get_el_amt_dict()
        cell_volume = self.structure.volume
        self._fraction_of_hydrogen_volume = ele_amt['H'] * 1.78 / cell_volume
        return self._fraction_of_hydrogen_volume

    @property
    def hydrogen_content(self):
        """
        calculate hydrogen content
        """
        ele_amt = self.structure.composition.get_el_amt_dict()
        self._hydrogen_content = ele_amt['H'] / self.structure.composition.num_atoms
        return self._hydrogen_content
    
    @property
    def density_of_hydrogen_per_volume(self):
        """
        The density of hydrogen atoms per unit volume.
            The volume of the total lattice minus the volume of the non-hydrogens is the volume occupied by the hydrogens
        """
        pass

    def describe(self, filename=None) -> pd.Series:
        """dict of key properties"""
        if self.hydride:
            return pd.Series({
                'formula': self.formula,
                'symmetry': self.spg_symmstry,
                'remain_H_ratio': self.remain_H_ratio,
                'fraction_of_hydrogen_volume': self.fraction_of_hydrogen_volume,
                'hydrogen_content': self.hydrogen_content,
                'check_dist': self.check_dist,
                'shr_num_avg': self.shr_num_avg,
                'cage_regularity_avg':  self.cage_regularity_avg,
                'h2h_network_regularity_avg': self.h2h_network_regularity_avg,
                'h2h_1nearest': self.h2h_1nearest[0][0],
                'h2h_2nearest': self.h2h_2nearest[0][0],
                'h2h_3nearest': self.h2h_3nearest[0][0],
                'nonh2nonh_1nearest' : self.nonh2nonh_1nearest[0][0],
                'nonh2nonh_2nearest' : self.nonh2nonh_2nearest[0][0],
                'hybrid': self.hybrid,
                'hydride': self.hydride,
                'filename': filename,
            })
        else:
            return pd.Series({
                'formula': self.formula,
                'symmetry': self.spg_symmstry,
                'remain_H_ratio': None,
                'fraction_of_hydrogen_volume': None,
                'hydrogen_content': None,
                'check_dist': None,
                'shr_num_avg': None,
                'cage_regularity_avg':  None,
                'h2h_network_regularity_avg': None,
                'h2h_1nearest':None,
                'h2h_2nearest': None,
                'h2h_3nearest': None,
                'nonh2nonh_1nearest' : None,
                'nonh2nonh_2nearest' : None,
                'hybrid': None,
                'hydride': None,
                'filename': filename,
            })
    


if __name__ == "__main__":
    from pymatgen.io.xyz import XYZ

    # f = Path('F:\\OneDrive - mails.jlu.edu.cn\\氢化物结构\\氢化物结构\\100_1.cif')
    # f = Path('./test/typical_hydrides/La2Y6H46_ContributedBy_Mayuan.vasp')
    # crystal = HydrideClathrate.from_file(f)
    # # print(crystal.hybrid_cage)
    # # cage = crystal.get_cage_mol(1)
    # # XYZ(cage).write_file("cage.xyz")
    # from pprint import pprint
    # pprint(crystal.h_network)
    # pprint(crystal.h2h_network_regularity)
    # pprint(crystal.h2h_network_regularity_avg)
    clathrate = HydrideClathrate.from_file("datasets/clathrate_hydrides/ThH10-Fm-3m-100GPa.cif")
    print(clathrate.hydrogen_content)
    
