import os
import sys
curPath = os.path.abspath(os.path.dirname(__file__))
rootPath = os.path.split(curPath)[0]
sys.path.append(rootPath)

import logging
from typing import List, Union

import numpy as np
from ase import Atoms

from .fingerprint import FingerPrint

logger = logging.getLogger("FP Mixins")


class FingerPrintMixin:
    """Base class to add fingerprint calculating method.

    All Optimizer should inherit this to set fingerprint to Atoms after
    optimized.
    """

    def __init__(self, fp_type: str, elems, threshold):
        self.fp_type = fp_type
        self.threshold = threshold
        self.fingerprint = FingerPrint(fp_type, elems, self.threshold)

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        return self.fingerprint.get_fp(atoms)

    def set_fp(self, atoms: Atoms):
        atoms.info['fp_type'] = self.fp_type
        atoms.info['fingerprint'] = self.get_fp(atoms)

    def get_similarity(self, fp1, fp2):
        dist, indexes, is_sim = self.fingerprint.get_similarity(fp1, fp2)
        return dist, indexes, is_sim


class UpdateBestMixin(FingerPrintMixin):
    """Base class to add pbest or gbest calculating method.

    All custom Structures should inherit this. Update best when new Atoms is
    added to the total structures
    """

    def __init__(self, critic: str, fp_type: str, elems, threshold: float):
        self.critic = critic
        self.threshold = threshold
        self._sort_func = self._get_sort_func()
        FingerPrintMixin.__init__(self, fp_type, elems, self.threshold)

    def update_best(self, structures: List[Atoms], nbest: int) -> List[Atoms]:
        structures = self.sorted_structures(structures)
        best_structures = list(self._sim_mat_filter(structures, nbest))
        return best_structures

    def sorted_structures(self, structures: List[Atoms]) -> List[Union[Atoms, None]]:
        '''
        首先使用 filter(self._fillout_nan, structures), 过滤掉焓值(enthalpy)为nan的Atoms结构, 将非nan的结构返回为List[Atoms]
        然后使用 sorted(List[Atoms], key=self._sort_func), 按照某一种标准对结构进行排序，这种标准有以下几种：
            _sort_enthalpy 按照焓值 对结构进行排序
            _sort_hardness 按照硬度 对结构进行排序
            ......
        '''
        structures = sorted(filter(self._fillout_nan, structures), key=self._sort_func)
        return structures

    def _sim_mat_filter(self, structures: List[Atoms], nbest: int):  # generator
        current_best = []
        for atoms in structures:
            if len(current_best) == 0 or not self.sim_condition(atoms, current_best):
                current_best.append(atoms)
                yield atoms
                if len(current_best) >= nbest:
                    return

    def _fillout_nan(self, atoms: Atoms) -> bool:
        '''
        过滤掉焓值(enthalpy)为nan的Atoms结构, 那么如何过滤呢？
        使用 self._sort_func(atoms) 对 atoms 结构进行过滤，根据选择的标准不同，
        _sort_func 函数会调用不同的形式：例如：
            _sort_enthalpy 按照焓值排除nan的形式
            _sort_hardness 按照硬度排除nan的形式
        '''
        x = self._sort_func(atoms)

        return not np.isnan(x)

    def _sort_enthalpy(self, atoms: Atoms) -> bool:
        return atoms.info.get('enthalpy', np.nan)

    def _sort_hardness(self, atoms: Atoms) -> bool:
        return -atoms.info('hardness', np.nan)

    def _get_sort_func(self):
        if self.critic.lower() == 'enthalpy':
            return self._sort_enthalpy
        elif self.critic.lower() == 'hardness':
            return self._sort_hardness
        elif self.critic.lower() == 'xrd':
            raise ValueError("Not Supported Yet!")
        else:
            raise ValueError(f"Unknown critic: {self.critic}!")

    def sim_condition(
        self, comparing_atoms: Atoms, compared_atoms_list: List[Atoms]
    ) -> bool:
        """Compare similarity of one Atoms to others

        Each atoms in `compared_atoms_list` is compared to `comparing_atoms`

        Args:
            comparing_atoms (Atoms): Atoms to compare.
            compared_atoms_list (List[Atoms]): Atoms list to be compared with.
            return_idx (bool, optional): Return the most similar index or not.
                Defaults to False.

        Returns:
            _type_: _description_
        """
        fp = comparing_atoms.info['fingerprint']
        compared_fps = np.array(
            [atoms.info['fingerprint'] for atoms in compared_atoms_list]
        )
        _, _, is_sim = self.fingerprint.get_similarity(fp, compared_fps)
        return is_sim
