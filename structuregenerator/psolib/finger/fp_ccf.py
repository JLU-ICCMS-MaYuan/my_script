import logging
from typing import List, Tuple, Union

import numpy as np
from ase import Atoms

from structuregenerator.psolib.finger.fp_base import BaseFingerPrint
from structuregenerator.psolib.lib_so.f90LegacyFingerPrint import Bondcrt

logger = logging.getLogger("CCF")


class CCF(BaseFingerPrint):
    def __init__(
        self,
        elems: List[str],
        threshold: float = 1e-4,
    ):
        self.elems = elems
        self.threshold = threshold

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        natom = len(atoms.numbers)  # total number of atoms
        unique, unique_indices, unique_counts = np.unique(
            atoms.numbers,
            return_index=True,
            return_counts=True,
        )
        ntype = len(unique)
        # matrix for f90wrapped function
        latmatrix = np.zeros((3, 3), dtype=np.float64, order='F')
        pos = np.zeros((3, natom), dtype=np.float64, order='F')
        type_num = np.zeros((ntype), dtype=np.int32, order='F')
        simf4 = np.zeros(
            (int((ntype * (ntype + 1)) / (2) * 400)), dtype=np.float64, order='F'
        )
        atomid = np.zeros((ntype + 1), dtype=np.int32, order='F')
        # set atomid, according to CryInitialGstruct.F90
        type_num = unique_counts[np.argsort(unique_indices)][:]
        atomid[0] = 1
        for i in range(ntype):
            atomid[i + 1] = atomid[i] + type_num[i]
        latmatrix = atoms.cell[:]
        pos = atoms.get_scaled_positions().T
        Bondcrt.storecfr(latmatrix, pos, atomid, simf4)
        ccf = simf4[:]
        return ccf

    def get_similarity(
        self,
        ccf1: Union[Atoms, np.ndarray],
        ccf2: Union[Atoms, List[Atoms], np.ndarray],
    ) -> Tuple[Union[np.ndarray, int, bool]]:

        if isinstance(ccf1, Atoms):
            ccf1 = self.get_fp(ccf1)
        if isinstance(ccf2, Atoms):
            ccf2 = [ccf2]
        if isinstance(ccf2, list):
            ccf2 = np.array([self.get_fp(atoms) for atoms in ccf2])
        if ccf2 is None:
            return 0, (0, 0), False
        dist = (ccf1 - ccf2)
        dist = np.linalg.norm(dist, axis=1)
        sim_idx = np.argmin(dist)
        nsim_idx = np.argmax(dist)
        return dist, (sim_idx, nsim_idx), dist[sim_idx] < self.threshold
