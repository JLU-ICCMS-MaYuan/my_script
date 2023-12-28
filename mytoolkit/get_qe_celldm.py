#!/usr/bin/env python3

import sys

from pathlib import Path

import numpy as np

# from pymatgen.core.structure import Structure
# a = s.lattice.a
# b = s.lattice.b
# c = s.lattice.c 
# alpha = s.lattice.alpha
# beta  = s.lattice.beta
# gamma = s.lattice.gamma
from ase.io import read 

try:
    filename = sys.argv[1]
except:
    filename = "POSCAR"

s = read(filename)
cell = s.cell.get_bravais_lattice()
a,b,c,alpha,beta,gamma = cell.cellpar()
crystal_family = cell.crystal_family
lattice_family = cell.lattice_system

print(a,b,c,alpha,beta,gamma)
print("crystal_family",crystal_family)
print("lattice_family",lattice_family)

celldm_1 = alat = a
celldm_2 = b/a
celldm_3 = c/a
celldm_4 = np.cos(np.pi/180*alpha) # cos(bc)
celldm_5 = np.cos(np.pi/180*beta) # cos(ac)
celldm_6 = np.cos(np.pi/180*gamma) # cos(ab)
print("celldm_1: {:<10.6f}\ncelldm_2: {:<10.6f}\ncelldm_3: {:<10.6f}\ncelldm_4: {:<10.6f}\ncelldm_5: {:<10.6f}\ncelldm_6: {:<10}\n".format(
    celldm_1,celldm_2,celldm_3,celldm_4,celldm_5,celldm_6
))

# celldm(2)-celldm(6)
# b,c,cosbc,cosac,cosab
# if lattice_family == 'cubic':
#     celldm_1 = alat = a
#     celldm_2 = 1 # b
#     celldm_3 = 1 # c
#     celldm_4 = 0 # cos(bc)
#     celldm_5 = 0 # cos(ac)
#     celldm_6 = 0 # cos(ab)
# elif lattice_family == 'hexagonal':
#     celldm_1 = alat = a
#     celldm_2 = b/a
#     celldm_3 = c/a
#     celldm_4 = np.cos(np.pi/180*alpha) # cos(bc)
#     celldm_5 = np.cos(np.pi/180*beta) # cos(ac)
#     celldm_6 = np.cos(np.pi/180*gamma) # cos(ab)
# elif lattice_family == 'hexagonal': # trigonal R c-axis
#     celldm_1 = alat = a
#     celldm_2 = b/a
#     celldm_3 = c/a
#     tx = np.sqrt(np.abs((1-c)/2))
#     ty = np.sqrt((np.abs(1-c)/6))
#     tz = np.sqrt((np.abs(1+2*c)/3))
#     v1 = a*np.array([tx, -ty, tz])
#     v2 = a*np.array([0, 2*ty, tz])
#     v3 = a*np.array([-tx, -ty, tz])
#     print(v1,v2,v3)
#     celldm_4 = celldm_5 = celldm_6 = np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
# elif lattice_family == 'rhombohedral' and crystal_family == 'hexagonal': # trigonal R c-axis
#     celldm_1 = alat = a
#     celldm_2 = b/a
#     celldm_3 = c/a
#     a = a/np.sqrt(3)
#     tx = np.sqrt(np.abs((1-c)/2))
#     ty = np.sqrt((np.abs(1-c)/6))
#     tz = np.sqrt((np.abs(1+2*c)/3))
#     u  = tz - 2*np.sqrt(2)*ty
#     v  = tz +   np.sqrt(2)*ty
#     v1 = a*np.array([u, v, v])
#     v2 = a*np.array([v, u, v])
#     v3 = a*np.array([v, v, u])
#     celldm_4 = celldm_5 = celldm_6 = np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))






