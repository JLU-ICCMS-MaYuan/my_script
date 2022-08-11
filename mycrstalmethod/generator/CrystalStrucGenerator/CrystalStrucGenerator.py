"""
generate 3D structure using method wrapped fortran code or pso
"""
import logging
from typing import List

import numpy as np
from ase import Atoms

from ...fingerprints.mixins import FingerPrintMixin
from ..base_generator import BaseStrucGenerator, is_bond_reasonable
from .f90sym3dgenerator import Spgcrylat, struct3dsym

logger = logging.getLogger("GENERATOR")


class CrystalStrucGenerator(BaseStrucGenerator, FingerPrintMixin):
    '''
    Crystal Structures Generators, which need to be initialized with config
    main function: generate(pbest, gbest, ispso)
    pbest and gbest must be delivered
    '''

    def __init__(self, config):
        self.config = config
        threshold = config['simthreshold']
        fingerprint: str = self.config['fingerprint']
        self.numberofspecies = len(self.config['nameofatoms'])
        elems = config['nameofatoms']
        FingerPrintMixin.__init__(self, fingerprint, elems, threshold)

    def generate_one(self, sid, pbest, gbest: List, current_atoms, fp_mats) -> Atoms:
        # if fp_mats is None:
        #     if self.fingerprint == 'ccf':
        #         fp_mats = np.zeros((2400, ))
        #     elif self.fingerprint == 'bcm':
        #         fp_mats = np.zeros((11, self.numberofspecies, self.numberofspecies))
        gen_method = self.choose_method(pbest, gbest, current_atoms)
        if gen_method == 'pso':
            gen_atoms = self.__pso_gen(pbest, gbest, current_atoms, fp_mats)
        elif gen_method == 'random':
            gen_atoms = self.__rand_gen(current_atoms, fp_mats)
        gen_atoms.info['sid'] = sid
        logger.info(f"Generating {sid} by {gen_method}, {gen_atoms.positions.shape}")
        return gen_atoms

    # pbest_*     :    [] or [Atoms]
    # pbest_list  :    [pbest_0, pbest_1, ...]
    # gbest       :    {'formula': [Atoms]}
    # -> list[Atoms]
    def generate_step(self, id_list, pbest_list, gbest, current_atoms_list, fp_mats):
        gen_structures_list = []
        for column, current_atoms in enumerate(current_atoms_list):
            gen_one = self.generate_one(
                id_list[column],
                pbest_list[column],
                gbest[current_atoms[0].get_chemical_formula('metal')],
                current_atoms,
                fp_mats,
            )
            gen_structures_list.append(gen_one)
        return gen_structures_list

    def choose_method(self, pbest, gbest: List, atoms):
        if len(pbest) == 0 and len(gbest) == 0:
            return 'random'
        elif len(pbest) * len(gbest) == 0 and self.config['runmode'] == 'sync_step':
            raise RuntimeError('either pbest or gbest exist but the other one is not')
        else:
            try:
                if atoms[1].info['next_evo_type'] == 'pso':
                    return 'pso'
                elif atoms[1].info['next_evo_type'] == 'random':
                    return 'random'
            #  in async mode, in first iteration, all strcut should be random
            #  but when gbest is not empty, atoms[1] will be initialized atoms
            #  which don't have'next evo type' attribution
            except KeyError:
                return 'random'

    def __pso_gen(
        self,
        pbest: List[Atoms],
        gbest: List[Atoms],
        current_atoms: List[Atoms],
        fp_mats,
    ):
        distance_of_ion = np.array(self.config['distanceofion'])
        init_atoms = current_atoms[0]
        opt_atoms = current_atoms[1]

        omega_max, omega_min = 0.9, 0.4
        iter_max = self.config['maxstep']
        step = current_atoms[0].info['sid'] // self.config['popsize']
        step = step + 1
        omega = omega_max - (omega_max - omega_min) * step / iter_max
        c1 = c2 = 2

        dists, indexes, _ = self.get_similarity(
            opt_atoms.info['fingerprint'],
            np.array([atoms.info['fingerprint'] for atoms in gbest]),
        )
        gbest = gbest[indexes[1]]
        pbest_pos = pbest[0].get_scaled_positions()
        # todo: decide which gbest should be used
        gbest_pos = gbest.get_scaled_positions()
        init_pos = init_atoms.get_scaled_positions()
        opt_volume = opt_atoms.get_volume()
        detla_pbest_init = pbest_pos - init_pos
        detla_gbest_init = gbest_pos - init_pos
        for _ in range(500):
            r1 = np.random.rand()
            r2 = np.random.rand()
            init_velocity = init_atoms.arrays['caly_velocity']
            velocity = (
                omega * init_velocity
                + c1 * r1 * detla_pbest_init
                + c2 * r2 * detla_gbest_init
            )
            gen_pso_pos = init_pos + velocity

            spacegroupid = np.random.uniform(1, 231, 1)
            lattice = np.ones((3, 3), dtype=np.float64, order='F')
            Spgcrylat().spggenlat(opt_volume, spacegroupid, lattice)
            # logger.warning("This lattice may not be physical.")

            new_atoms = Atoms(
                symbols=opt_atoms.symbols,
                scaled_positions=gen_pso_pos,
                cell=lattice,
                pbc=True,
            )
            new_atoms.arrays['caly_velocity'] = velocity
            fp = self.get_fp(new_atoms)
            name_of_atoms = self.config['nameofatoms']
            if (
                is_bond_reasonable(new_atoms, name_of_atoms, distance_of_ion)
                and not self.get_similarity(fp, fp_mats)[2]  # not is_sim
            ):
                break
        else:
            logger.warning("This structure may not be physical.")

        self.set_fp(new_atoms)
        return new_atoms

    def __rand_gen(self, current_atoms, fp_mats):
        symbols = current_atoms[0].symbols
        total_numb_atoms = len(current_atoms[0])

        name_of_atoms = self.config['nameofatoms']
        number_of_atoms = np.array(
            [symbols.formula.count()[name] for name in name_of_atoms],
            dtype=np.int32,
            order='F',
        )

        volume_per_atom = np.array(self.config['volumeperatom'])
        volume = sum(volume_per_atom * number_of_atoms)

        distance_of_ion = self.config['distanceofion']
        distance_of_ion = np.array(distance_of_ion, dtype=np.float64, order='F')
        # temp
        ocmatrix = np.ones((27, len(name_of_atoms)), dtype=np.int32, order='F')

        # out
        lattice = np.ones((3, 3), dtype=np.float64, order='F')
        position = np.ones((3, total_numb_atoms), dtype=np.float64, order='F')
        # velocity = np.ones((3, total_numb_atoms), dtype=np.float64, order='F')

        for i in range(100):
            ocmatrix = np.ones((27, len(name_of_atoms)), dtype=np.int32, order='F')

            lattice = np.ones((3, 3), dtype=np.float64, order='F')
            position = np.ones((3, total_numb_atoms), dtype=np.float64, order='F')

            # velocity = np.ones((3, total_numb_atoms), dtype=np.float64, order='F')
            spacegroup_id, lgen = struct3dsym(
                len(name_of_atoms),
                number_of_atoms,
                total_numb_atoms,
                volume,
                distance_of_ion,
                self.config['spespacegroup'][0],
                self.config['spespacegroup'][1],
                ocmatrix,
                lattice,
                position,
            )
            # gfortran version
            position = position.T
            # velocity = velocity.T  # is zeros
            velocity = np.random.uniform(0, 1, (total_numb_atoms, 3))
            if lgen == 1:
                random_atoms = Atoms(
                    symbols=symbols,
                    scaled_positions=position,
                    cell=lattice,
                    pbc=True,
                )
                fp = self.get_fp(random_atoms)
                if (
                    is_bond_reasonable(random_atoms, name_of_atoms, distance_of_ion)
                    and not self.get_similarity(fp, fp_mats)[2]  # not is_sim
                ):
                    break
        else:
            logger.warning("Random structure may not be physical")
        random_atoms.arrays['caly_velocity'] = velocity
        self.set_fp(random_atoms)
        return random_atoms
