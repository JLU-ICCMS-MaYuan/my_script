"""Sort Atoms element sequence.

pervoski = Atoms(...)
pervoski.get_chemical_symbols() -> ['Ca', 'Ti', 'O', 'O', 'O']
pervoski = sort_atoms(pervoski, ['Ti', 'O', 'Ca']
pervoski.get_chemical_symbols() -> ['Ti', 'O', 'O', 'O', 'Ca']
"""

from itertools import chain
from typing import List, Union

from ase import Atoms
from ase.formula import Formula


def sort_atoms(atoms: Union[Atoms, str, Formula], elems: List) -> Atoms:
    if isinstance(atoms, Atoms):
        elems_d = {elem: [] for elem in elems}
        for idx, atom in enumerate(atoms):
            elems_d[atom.symbol].append(idx)
        indices = list(chain.from_iterable(elems_d.values()))
        return atoms[indices]
    elif isinstance(atoms, str):
        f_count = Formula(atoms).count()
        return ''.join(elem + str(f_count[elem]) for elem in elems)
    elif isinstance(atoms, Formula):
        f_count = atoms.count()
        return Formula(''.join(elem + str(f_count[elem]) for elem in elems))
    else:
        raise ValueError('Only accept {Atoms|str|Formula}')
