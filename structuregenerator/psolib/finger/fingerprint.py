import logging

import numpy as np
from ase import Atoms

from . import fp_bcm
from . import fp_ccf

logger = logging.getLogger("Fingerprint")


class FingerPrint:
    def __init__(self, fingerprint: str, elems, threshold):
        fp_types = ['bcm', 'ccf', 'dummy']
        self.threshold = threshold
        if fingerprint not in fp_types:
            raise ValueError(
                f"Unknown fingerprint: {fingerprint}! Should be {fp_types}"
            )
        elif fingerprint == 'bcm':  # bond charactor matrix
            py = False
            if py:
                self.calculator = fp_bcm.BondCharMatrix(elems, threshold=self.threshold)
            else:
                self.calculator = fp_bcm.BondCharMatrix_F90(
                    elems, threshold=self.threshold
                )
        elif fingerprint == 'ccf':
            self.calculator = fp_ccf.CCF(elems, threshold=self.threshold)
        elif fingerprint == 'dummy':
            pass
        else:
            raise RuntimeError("FingerPrint BUG! Please Report or Fix")

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        return self.calculator.get_fp(atoms)

    def get_similarity(self, fp1, fp2):
        dist, indexes, is_sim = self.calculator.get_similarity(fp1, fp2)
        return dist, indexes, is_sim
