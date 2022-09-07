#!/usr/bin/env python
import numpy as np


def uniq_atom_names(data):
    """
    Make the atom names uniq. For example
    ['O', 'H', 'O', 'H', 'O'] -> ['O', 'H']
    Parameters
    ----------
    data : dict
        data dict of `System`, `LabeledSystem`
    """
    unames = []
    uidxmap = []
    for idx, ii in enumerate(data['atom_names']):
        if ii not in unames:
            unames.append(ii)
        uidxmap.append(unames.index(ii))
    data['atom_names'] = unames
    tmp_type = list(data['atom_types']).copy()
    data['atom_types'] = np.array([uidxmap[jj] for jj in tmp_type], dtype=int)
    data['atom_numbs'] = [
        sum(ii == data['atom_types']) for ii in range(len(data['atom_names']))
    ]
    return data
