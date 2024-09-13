import abc
from typing import Tuple, Union

import numpy as np
from ase import Atoms


class BaseFingerPrint(abc.ABC):
    @abc.abstractmethod
    def get_fp(self, atoms: Atoms) -> np.ndarray:
        pass

    @abc.abstractmethod
    def get_similarity(self, fp1, fp2) -> Tuple[Union[np.ndarray, Tuple[int], bool]]:
        """Get similarity distance between two structures or fp1 to all fp2.

        If fp1 is a Atoms, fp1 should be turned into ndarray. If fp2 is one
            Atoms or a list of Atoms, fp2 should be turned into ndarray with
            first dim of number of Atoms.

        Return distance of fp1 to fp2, and tuple of the most similar and not
            similar indexes in fp2, and if these two are similar or not
            (smaller than threshold).

        Args:
            fp1 : Atoms or fingerprint.
            fp2 : Atoms or list of Atoms or fingerprints.

        Returns:
            Tuple[Union[np.ndarray, Tuple[int], bool]]: dists, (min, max), is_sim
        """
        pass
