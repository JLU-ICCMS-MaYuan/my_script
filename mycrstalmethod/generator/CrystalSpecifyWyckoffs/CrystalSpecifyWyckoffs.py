import logging
import random
import re
from collections import Counter, defaultdict
from itertools import chain, combinations, product

import numpy as np
from ase import Atoms
from ase.formula import Formula

try:
    from pyxtal import pyxtal
except (ImportError, ModuleNotFoundError):
    no_pyxtal = True
else:
    no_pyxtal = False

from ...utils.sort_atoms import sort_atoms
from ..base_generator import is_bond_reasonable
from ...fingerprints.mixins import FingerPrintMixin
from ..CrystalStrucGenerator.CrystalStrucGenerator import CrystalStrucGenerator

logger = logging.getLogger("GENERATOR")


class CrystalSpecifyWyckoffs(CrystalStrucGenerator):
    def __init__(self, config):
        if no_pyxtal:
            raise ModuleNotFoundError(
                "No package 'pyxtal', please try `pip install pyxtal`"
            )
        self.config = config
        threshold = config['simthreshold']
        fingerprint: str = self.config['fingerprint']
        self.numberofspecies = len(self.config['nameofatoms'])
        self.nameofatoms = self.config['nameofatoms']
        elems = config['nameofatoms']
        FingerPrintMixin.__init__(self, fingerprint, elems, threshold)
        self._group = self.get_group(
            config['optionalsites'], config['sitesoccupiedrange']
        )

    # generate_one  same as CrystalStrucGenerator.generate_one
    # generate_step same as CrystalStrucGenerator.generate_step
    # choose_method same as CrystalStrucGenerator.choose_method
    # _pso_gen      same as CrystalStrucGenerator._pso_gen

    def __rand_gen(self, current_atoms, fp_mats) -> Atoms:
        '''
        create struct by choosing wyckoff positions
        '''
        spespacegroup = self.config['spespacegroup'][0]
        distance_of_ion = self.config['distanceofion']
        name_of_atoms = self.config['nameofatoms']
        symbols = current_atoms[0].symbols.get_chemical_formula('metal')

        species_amounts_sites = self._group[symbols]
        spe_amo_site = random.choice(species_amounts_sites)

        amounts = spe_amo_site[0]
        wyck = spe_amo_site[1]
        struc = pyxtal()
        for i in range(100):
            try:
                struc.from_random(
                    3,
                    spespacegroup,
                    name_of_atoms,
                    amounts,
                    factor=2.0,
                    sites=wyck,
                    # lattice=my_lat
                )
                total_num_atoms = sum(amounts)

                _random_atoms = struc.to_ase()
                random_atoms = sort_atoms(_random_atoms, name_of_atoms)

                fp = self.get_fp(random_atoms)
                if (
                    is_bond_reasonable(random_atoms, name_of_atoms, distance_of_ion)
                    and not self.get_similarity(fp, fp_mats)[2]  # not is_sim
                ):
                    break
            except Exception as e:
                logger.debug(f"except {e}")
                continue
        else:
            logger.warning("this structure may not be physical")

        velocity = np.random.uniform(0, 1, (total_num_atoms, 3))
        random_atoms.arrays['caly_velocity'] = velocity
        self.set_fp(random_atoms)
        return random_atoms

    def get_H(self, H_occupied_wps, h_lower, h_upper):
        if h_upper < h_lower:
            h_lower, h_upper = h_upper, h_lower
        hydrogen_wps = []
        for i in range(h_lower, h_upper + 1):
            for e in product(H_occupied_wps, repeat=i):
                e = sorted(list(e))
                res = Counter(e)
                if (np.array(list(res.values())) < 2).all():
                    if not hydrogen_wps:
                        hydrogen_wps.append(e)
                    else:
                        for h_wp in hydrogen_wps:
                            h_wp = sorted(list(h_wp))
                            if h_wp == e:
                                break
                        else:
                            hydrogen_wps.append(e)

        hydrogen_atoms_number = []
        for wps_strings in hydrogen_wps:
            atoms = []
            for wp_multi in wps_strings:
                atoms.extend(re.findall(r"\d+", wp_multi))
            atoms = list(map(int, atoms))
            hydrogen_atoms_number.append(sum(atoms))

        return hydrogen_wps, hydrogen_atoms_number

    def get_M(self, X_occupied_wps, X_lower, X_upper):
        if X_upper < X_lower:
            X_lower, X_upper = X_upper, X_lower
        X_wps = []
        for i in range(X_lower, X_upper + 1):
            for e in combinations(X_occupied_wps, i):
                e = sorted(list(e))
                X_wps.append(e)
        X_atoms_number = []
        for wps_strings in X_wps:
            atoms = []
            for wp_multi in wps_strings:
                atoms.extend(re.findall(r"\d+", wp_multi))
            atoms = list(map(int, atoms))
            X_atoms_number.append(sum(atoms))
        return X_wps, X_atoms_number

    def get_X(self, X, X_occupied_wps, X_lower, X_upper):
        if X == 'H':
            return self.get_H(X_occupied_wps, X_lower, X_upper)
        else:
            return self.get_M(X_occupied_wps, X_lower, X_upper)

    def get_group(self, optionalsites, sitesoccupiedrange) -> dict:
        """Function

        Args:
            optional_sites (List[List[str]]):
                optional wyckoff positions of each elements,
                e.g. [['4a', "6b", "8c"],['12g', '12h', "16g"]],
                     ['12m', '16g'].
            sitesoccupiedrange (List[List[int]]):
                repeat range of each wyckoff position,
                e.g. [[1, 3], [2, 6], [2, 6]].

        Returns:
            dict: {'formula0': [amount, wyck], ...}
        """
        elems = self.nameofatoms  # ['H', 'Ca', 'C']
        _group = defaultdict(list)
        wyck_list, nelems_list = [], []
        for elem, opt_sites, sites_range in zip(
            elems, optionalsites, sitesoccupiedrange
        ):
            wyck, nelems = self.get_X(elem, opt_sites, *sites_range)
            wyck_list.append(
                wyck
            )  # [A, B, C]: [[['4b'], ['4d', '4f']], [['4b'], ['4d', '4f'], ...]
            nelems_list.append(
                nelems
            )  # [A, B, C]: [[4, 8], [4, 8], ...]
        nelems_comb = product(*nelems_list)  # generator
        wyck_comb = product(*wyck_list)
        for nelems, wyck in zip(nelems_comb, wyck_comb):
            formula = ''.join(map(str, chain.from_iterable(zip(elems, nelems))))
            formula = Formula(formula).format("metal")
            _group[formula].append([nelems, wyck])
        _group = dict(_group)
        return _group
