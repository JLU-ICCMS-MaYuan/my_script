#!/usr/bin/env python

import logging
import itertools

from ase import Atoms

from structuregenerator.psolib.utils.sort_atoms import sort_atoms

logger = logging.getLogger("convert")


def dict2Atoms(data, references_atom_names):
    atom_names = data['atom_names']
    atom_numbs = data['atom_numbs']
    cells = data['cells'][-1]
    coords = data['coords'][-1]
    _symbols = list(itertools.chain.from_iterable(zip(atom_names, atom_numbs)))
    symbols = ''.join(list(map(str, _symbols)))
    _atoms = Atoms(
        symbols=symbols,
        positions=coords,
        cell=cells,
    )
    atoms = sort_atoms(_atoms, references_atom_names)
    return atoms
