#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

import spglib
import cellconstructor as CC

def get_symmetry(number_of_population, irr_idx, precision):
    space_groups = []
    angles = []
    paras = []
    for dyn_pop_idx in range(number_of_population):
        # Load the dynamical matrix
        dyn = CC.Phonons.Phonons(f"dyn_pop{dyn_pop_idx}_", nqirr = irr_idx)
        # Get the structure
        structure = dyn.structure.get_ase_atoms()
        unit_cell = structure.unit_cell
        # Get the symmetry information
        spgroup = spglib.get_spacegroup(structure, 0.05)
        space_groups.append(spgroup)

        para = np.sqrt(unit_cell[0,0]*unit_cell[0,0]+unit_cell[0,1]*unit_cell[0,1]+unit_cell[0,2]*unit_cell[0,2])
        paras.append(para)
        np.savetxt("para.dat", paras)

        angle = np.zeros(3)
        for i in range(3):
            nexti = (i+1)%3
            otheri = (i+2)%3
            angle[otheri] = np.arccos( np.dot(unit_cell[i,:], unit_cell[nexti,:]) /
                (np.sqrt(np.dot(unit_cell[i,:], unit_cell[i,:])) *
                np.sqrt(np.dot(unit_cell[nexti,:], unit_cell[nexti,:])))) * 180 / np.pi
        angles.append(angle)
        np.savetxt("angle.dat", angles)
        
        # We can save them in the output at each minimization step
        f = open("space_group.dat", "w")
        f.writelines(["{}) {}\n".format(i+1, x) for i,x in enumerate(space_groups)])
        f.close()

if __name__ == '__main__':
    parser = ArgumentParser(description="Convert QE dynamical matrix file to POSCAR format.\n You can use this scrips by:   qedyn2struct.py -i ../V3_Hessian.dyn1 -o qe -kd 8 8 8 -ks 4 4 4 -pp /work/home/mayuan/workplace/5.calypso/35.Ce-Sc-H/4.detailed-compute/1.CeSc2H24/200GPa/0.prepare-init-dyn/pp -pn Ce Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf Sc Sc.pbe-spn-kjpaw_psl.1.0.0.UPF H H.pbe-kjpaw_psl.1.0.0.UPF")
    parser.add_argument('-n', '--number-of-population', type=int, default=0, help="input number of population")
    parser.add_argument('-i', '--irr_idx', type=int, required=True, help="input irr_idx")
    parser.add_argument('-p', '--precision', type=float, default=0.01, help="Input precision for find symmetry")
    
    args = parser.parse_args()
    number_of_population = args.number_of_population
    irr_idx = args.irr_idx
    precision = args.precision
    
    get_symmetry(number_of_population, irr_idx, precision)

