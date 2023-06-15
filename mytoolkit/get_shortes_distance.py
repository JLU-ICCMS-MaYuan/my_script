#!/usr/bin/env python

from ase.io import read

import numpy as np

atom = read("POSCAR")
d = atom.get_all_distances()
np.fill_diagonal(d, 10000)
short_d = np.min(d)
print(short_d)