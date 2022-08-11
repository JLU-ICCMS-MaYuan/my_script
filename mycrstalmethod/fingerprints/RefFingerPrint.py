import numpy as np
from ase import Atoms


class FingerPrint:
    def __init__(self, fingerprint: str):
        fp_types = ['bcm', 'ccf', 'dummy']
        if fingerprint not in fp_types:
            raise ValueError(
                f"Unknown fingerprint: {fingerprint}! Should be {fp_types}"
            )
        elif fingerprint == 'bcm':  # bond character matrix
            self.calc_fp = self._calc_bcm
        elif fingerprint == 'ccf':
            self.calc_fp = self._calc_ccf
        elif fingerprint == 'dummy':
            self.calc_fp = self._calc_dummy
        else:
            raise RuntimeError("FingerPrint BUG! Please Report or Fix")

    def get_fp(self, atoms: Atoms) -> np.ndarray:
        return self.calc_fp(atoms)

    # Bond Character Matrix
    def _calc_bcm(self, atoms: Atoms) -> np.ndarray:
        from FingerPrintsGenerator.FingerPrints import Bondcharmatrix

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
        simf1 = np.zeros((11, ntype, ntype), dtype=np.float64, order='F')
        latmatrix = atoms.cell[:, :]
        pos = atoms.get_scaled_positions().T
        type_num = unique_counts[np.argsort(unique_indices)][:]
        Bondcharmatrix.bondcharmatrix_calbcm(
            1, ntype, natom, type_num, latmatrix, pos, simf1
        )
        bcm = simf1[:, :, :]
        return bcm

    def _calc_ccf(self, atoms: Atoms) -> np.ndarray:
        from FingerPrintsGenerator.FingerPrints import Bondcrt

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

    def _calc_dummy(self, atoms: Atoms) -> np.ndarray:
        return 0.0


if __name__ == '__main__':
    """
    ! Calculate the distance between the optimized structure of AFM and FM magnetic type
    !
    """
    import argparse
    import numpy as np
    from matplotlib import pyplot as plt
    import ase.io
    import os

    plt.style.use(['science', 'ieee', 'no-latex'])

    # Create the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--WorkDir", type=str, required=True)
    args = parser.parse_args()

    myFingerprint = FingerPrint(fingerprint="ccf")
    StructDiff = []
    AFMccf = []
    FMccf = []
    for name in ["CoCo", "CrCr", "FeFe", "MoMo", "TiTi", "VV"]:
        atoms = ase.io.read(os.path.join(args.WorkDir, name, "AFM", "CONTCAR"))
        AFMccf.append(myFingerprint.get_fp(atoms))
        atoms = ase.io.read(os.path.join(args.WorkDir, name, "FM", "CONTCAR"))
        FMccf.append(myFingerprint.get_fp(atoms))
        # StructDiff.append(np.linalg.norm(FMccf - AFMccf))
    for i in range(len(AFMccf)):
        for j in range(len(AFMccf)):
            StructDiff.append(np.linalg.norm(FMccf[i] - AFMccf[j]))
    print(StructDiff)

    from matplotlib import pyplot as plt

    plt.style.use(['science', 'ieee', 'no-latex'])

    plt.figure()
    # plt.xticks(np.arange(6), ["CoCo", "CrCr", "FeFe", "MoMo", "TiTi", "VV"])
    plt.plot(StructDiff, "bs-")
    # plt.ylabel(r"$Distance$")
    # plt.ylim((0, 100))
    plt.savefig(os.path.join(args.WorkDir, "DistanceTemp.png"))
