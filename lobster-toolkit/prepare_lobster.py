#!/usr/bin/env python3

from ase.io.pov import get_bondpairs

bondpairs_raw = get_bondpairs(atoms)
bondpairs = []
for bp in bondpairs_raw:
    tmp = list(bp[:2])
    tmp.sort()
    if tmp not in bondpairs and tmp[0] != tmp[1]:
        bondpairs.append(tmp)