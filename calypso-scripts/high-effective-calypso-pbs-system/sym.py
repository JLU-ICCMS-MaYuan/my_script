#!/data/home/mym/soft/anaconda/ana/bin/python
import os
import glob
import time
import math
import numpy as np
import libpyspg as spg
def structure():
 fline = []
 f = open('POSCAR')
 for line in f:
    fline.append(line)
 lat = []
 for i in range(2, 5):
     lat.append(map(float, fline[i].split()))
 lat = np.array(lat, float)
# print lat
 typt = map(int, fline[6].split())
 natom=np.array(sum(typt),int)
# print typt
 pos = []
 for item in fline[8:8+natom]:
     pos.append(map(float, item.split()[:3]))
 pos = np.array(pos, float)
 f.close
# print pos
 cell=[]
 cell.append(lat)
 cell.append(pos)
 cell.append(typt)
 cell.append(natom)
 spgdata = findsym(cell, 0.5)
 print spgdata
def findsym(cell, prec):
    '''
    cell[ lattice,
          positions,
          typt,
          num_atoms
    ]
    '''

    aprec = -1

    sl = cell[0].T.copy()
    sp = cell[1].copy()
    stypt = cell[2][:]
    num_atoms = cell[3]

    snumbers = []
    for i in range(len(stypt)):
        snumbers += [i+1] * stypt[i]
    snumbers = np.array(snumbers, int)

    l = sl[:]
    p = sp[:]
    typt = stypt[:]
    numbers = snumbers[:]

    (num_spg, symbol_spg) = spg.spacegroup(l, p, numbers, prec, aprec)
    return num_spg
if __name__ == '__main__':
 structure()
