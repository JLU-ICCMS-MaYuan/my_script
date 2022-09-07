import logging
from itertools import product
from typing import List, Tuple, Union

import numpy as np
from ase import Atoms
from scipy.special import sph_harm

from psolib.utils.sort_atoms import sort_atoms
from psolib.finger.fp_base import BaseFingerPrint
from psolib.lib_so.f90LegacyFingerPrint import Bondcharmatrix

logger = logging.getLogger("BCM")

def cart2sph(cart: np.ndarray) -> np.ndarray:
    """Cartesian to Spherical

    (*, (x, y, z)) -> (*, (r, theta, phi))  unit in rad

    Transform on the last dimensions.

    Input an array with cartesian coordinates in last dim, shape of (*, 3),
    transform to spherical coordinates. Return shape of (*, 3), last dim of
    (radius theta phi). Keep the other dimensions.

    Args:
        cart (np.ndarray): Cartesian coordinates.

    Returns:
        np.ndarray: Spherical coordinates.
    """
    cart = np.moveaxis(cart, -1, 0)
    xy = cart[0] ** 2 + cart[1] ** 2
    r = np.sqrt(xy + cart[2] ** 2)
    theta = np.arctan2(np.sqrt(xy), cart[2])
    phi = np.arctan2(cart[1], cart[0])
    sph = np.stack([r, theta, phi], 0)
    sph = np.moveaxis(sph, 0, -1)
    return sph


class BondCharMatrix(BaseFingerPrint):
    """Calculate Bond Charactor Matrix. <doi:10.1016/j.cpc.2012.05.008>

    Q_(lm)^(AB) = 1/N_AB * sum_(i,j)[exp(-alpha * (r_ij - b_AB)) *
    Y_lm(theta_ij, phi_ij)]

    Q_l^(AB) = sqrt{4pi/(2l + 1) + sum_m[|Q_(lm)^(AB)|^2]}

    Split pair distances into block, ordered by elements.
    [A, B, C] -> [AA, AB, AC, BA, BB, BC, CA, CB, CC]  (itertools.product)

            A1  A2  B1  B2
        A1    A-A  |  A-B
        A2 ________|______
        B1    B-A  |  B-B
        B2         |

    Args:
        elems (List[str], optional): List of elements, to sort atoms in order.
            If not given, will be set to the sequence found in calculated
            structure, e.g. ['Li', 'Na'].
        rcut (float, optional): Cutoff radius of bond. Defaults to 3.5.
        alpha (float, optional): Radius coefficient. Defaults to 0.
        max_l (int, optional): Maximum angular momentum number. Defaults to 11.
        threshold (float, optional): Similarity threshold. Defaults to 1e-4.
    """

    def __init__(
        self,
        elems: List[str],
        rcut: float = 3.5,
        alpha: float = 0,
        max_l: int = 11,
        threshold: float = 1e-4,
    ):
        self.elems = elems
        self.rcut = rcut
        self.alpha = alpha
        self.max_l = max_l
        self.threshold = threshold

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        """Get BCM.

        Element sequence MUST be given.

        Args:
            atoms (Atoms): Input structure

        Returns:
            np.ndarray: Bond charactor matrix, (Max_L, NumberOfElem, NumberOfElem).
        """
        # sort atoms and get number of each element
        elems = self.elems
        atoms = sort_atoms(atoms, elems)
        formula_count = atoms.symbols.formula.count()  # a dict {'A': nA, 'B': nB}
        nelems = list(formula_count.values())  # number of each element
        # get distance vector (natoms, natoms, 3)
        all_dists = atoms.get_all_distances(mic=True, vector=True)
        sph_all_dists = cart2sph(all_dists)  # (natoms, natoms, 3), (r, theta, phi)
        eps = 1e-1  # minimum bond length to be considered
        Ql = np.zeros([self.max_l, len(elems), len(elems)])
        for angular_l in range(0, self.max_l):  # Angular momentum number
            if (
                angular_l % 2 == 0
            ):  # only even-l are invariant with respect to direction
                for m in range(-angular_l, angular_l + 1):  # Magnetic momentum number
                    Qlm = np.zeros(
                        [2 * angular_l + 1, len(elems), len(elems)],
                        dtype=np.complex64,
                    )
                    for i, j in product(
                        range(len(elems)), repeat=2
                    ):  # split into block
                        block_sph = sph_all_dists[
                            sum(nelems[:i]) : sum(nelems[: i + 1]),
                            sum(nelems[:j]) : sum(nelems[: j + 1]),
                        ]  # each A-B block, [n, n, 3]
                        block_sph = block_sph.reshape(-1, 3)  # flatten, (N, 3)
                        # filter bond length between [eps, rcut]
                        bonded_sph = [di for di in block_sph if eps < di[0] < self.rcut]
                        if len(bonded_sph) == 0:  # no bonds between A-B
                            Qlm[m, i, j] = 0
                        else:
                            bonded_sph = np.atleast_2d(bonded_sph)  # check shape
                            bAB = np.min(bonded_sph[:, 0])  # min bond length of A-B
                            Y = sph_harm(
                                m, angular_l, bonded_sph[:, 1], bonded_sph[:, 2]
                            )
                            exp_part = np.exp(-self.alpha * (bonded_sph[:, 0] - bAB))
                            Qlm[m, i, j] = np.mean(exp_part * Y)
                q_power = np.sum(np.power(np.abs(np.real(Qlm)), 2), axis=0)
                nom_coef = 4 * np.pi / (2 * angular_l + 1)
                Ql[angular_l] = np.sqrt(nom_coef * q_power)
        return Ql

    def get_similarity(
        self,
        bcm1: Union[Atoms, np.ndarray],
        bcm2: Union[Atoms, List[Atoms], np.ndarray],
    ) -> Tuple[Union[np.ndarray, int, bool]]:
        """Get similarity distance between atoms1 and atoms2

        Atoms1 and Atoms2 must have the same element order. bcm1 is shaped as
            [k, k], bcm2 is reshaped to [N, k, k].

        Return distances of shape [N], the most similar index (in N), and is
            the smallest distance smaller than threshold or not.

        Args:
            bcm1 (Union[Atoms, np.ndarray]): Structure1 or bcm1.
            bcm2 (Union[Atoms, List[Atoms], np.ndarray]): Structure2 or bcm2.

        Returns:
            Tuple[Union[np.ndarray, int, bool]]: Distance, bool of is similar
                or not, and most similar index.
        """
        if isinstance(bcm1, Atoms):
            bcm1 = self.get_fp(bcm1)
        if isinstance(bcm2, Atoms):
            bcm2 = [bcm2]
        if isinstance(bcm2, list):
            bcm2 = np.array([self.get_fp(atoms) for atoms in bcm2])
        if bcm2 is None:
            return 0, (0, 0), False
        if bcm2.ndim != 4:
            bcm2 = np.reshape(bcm2, (-1, self.max_l, len(self.elems), len(self.elems)))
        dist = (bcm1 - bcm2).reshape(bcm2.shape[0], -1)
        dist = np.linalg.norm(dist, axis=1)
        sim_idx = np.argmin(dist)
        nsim_idx = np.argmax(dist)
        return dist, (sim_idx, nsim_idx), dist[sim_idx] < self.threshold


class BondCharMatrix_F90(BaseFingerPrint):

    def __init__(
        self,
        elems: List[str],
        max_l: int = 11,
        threshold: float = 1e-4,
    ):
        self.elems = elems
        self.max_l = max_l
        self.threshold = threshold

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        natom = len(atoms.numbers)
        unique, unique_indices, unique_counts = np.unique(
            atoms.numbers, return_index=True, return_counts=True
        )
        ntype = len(unique)
        # matrix for f90wrapped function
        latmatrix = np.zeros((3, 3), dtype=np.float64, order='F')
        pos = np.zeros((3, natom), dtype=np.float64, order='F')
        type_num = np.zeros((ntype), dtype=np.int32, order='F')
        simf1 = np.zeros((11, ntype, ntype), dtype=np.float64, order='F')
        latmatrix = atoms.cell[:, :]
        pos = atoms.get_scaled_positions().T
        type_num = unique_counts[np.argsort(unique_indices)][:]
        Bondcharmatrix.bondcharmatrix_calbcm(1, ntype, natom, type_num, latmatrix, pos, simf1)
        bcm = simf1[:, :, :]
        return bcm

    def get_similarity(
        self,
        bcm1: Union[Atoms, np.ndarray],
        bcm2: Union[Atoms, List[Atoms], np.ndarray],
    ) -> Tuple[Union[np.ndarray, int, bool]]:

        if isinstance(bcm1, Atoms):
            bcm1 = self.get_fp(bcm1)
        if isinstance(bcm2, Atoms):
            bcm2 = [bcm2]
        if isinstance(bcm2, list):
            bcm2 = np.array([self.get_fp(atoms) for atoms in bcm2])
        if bcm2 is None:
            return 0, (0, 0), False
        if bcm2.ndim != 4:
            bcm2 = np.reshape(bcm2, (-1, self.max_l, len(self.elems), len(self.elems)))
        dist = (bcm1 - bcm2).reshape(bcm2.shape[0], -1)
        dist = np.linalg.norm(dist, axis=1)
        sim_idx = np.argmin(dist)
        nsim_idx = np.argmax(dist)
        return dist, (sim_idx, nsim_idx), dist[sim_idx] < self.threshold
