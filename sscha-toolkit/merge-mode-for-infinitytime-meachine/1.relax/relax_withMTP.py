from __future__ import print_function
import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.calculators

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Cluster
import spglib

from numpy import *
import numpy as np

import os

def convert_dat2cfg(population_id):
    # Generates the CFG files used by the MTP 
    cfg_file = open("MLIP_"+str(population_id) + ".cfg", "w")
    directory = f"ensemble_dat"
    FORCE_FILES = True

    energy_file = open(directory + "/energies_supercell_population" + str(population_id) + ".dat", "r")
    energy_lines = [l.strip() for l in energy_file.readlines()]
    energy_file.close()

    print(len(energy_lines),energy_lines)

    for j in range(1,len(energy_lines)+1):
        structure_file = open(directory + "/scf_population" + str(population_id) + "_" + str(j) + ".dat", "r")
        lines = [l.strip() for l in structure_file.readlines()]
        structure_file.close()

        try:
            force_file = open(directory + "/forces_population" + str(population_id) + "_" + str(j) + ".dat", "r")
            force_lines = [l.strip() for l in force_file.readlines()]
            force_file.close()
        except: FORCE_FILES = False

        try:
            pressure_file = open(directory + "/pressures_population" + str(population_id) + "_" + str(j) + ".dat", "r")
            pressure_lines = [l.strip() for l in pressure_file.readlines()]
            pressure_file.close()
        except:
            pressure_lines = [
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            ]

        SIZE = len(lines)-6
        cell = np.zeros((3, 3))
        atoms = np.zeros((SIZE,3))
        forces = np.zeros((SIZE,3))
        atm_type = [None] * SIZE
        pressures = np.zeros((3, 3))

        for i in range(0,3):
            cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]

        for i in range(0,SIZE):
            atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
            atm_type[i] =  lines[i+6].split()[0]
            if FORCE_FILES == True:
                forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                forces[i,:] = forces[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__


        for i in atm_type:
            try: self.atom_type_table[i]
            except: self.atom_type_table[i] = len(self.atom_type_table)


        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))

        for i in range(0,3):

            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]
            pressures[i, :] = pressures[i, :] * self.__Ry_to_eV__ / self.__Bohr_to_A__ / self.__Bohr_to_A__ / self.__Bohr_to_A__* Volume

        cfg_file.write("BEGIN_CFG\n")
        cfg_file.write(" Size\n")
        cfg_file.write("{: >5}".format(SIZE) +"\n")
        cfg_file.write(" Supercell\n")
        for row in cell:
            cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
        cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for i in range(0,len(atoms)):
            if FORCE_FILES == True:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(self.atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
            else:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(self.atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*[0,0,0]))
        cfg_file.write(" Energy\n")
        cfg_file.write("     {: >13f}".format(float(energy_lines[j-1])*self.__Ry_to_eV__) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(self.atom_type_table) + "\n")
        cfg_file.write(" Feature conf_number " + str(j) + "\n")
        cfg_file.write(" Feature population " + str(population_id) + "\n")
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()

def convert_cfg2dat(population_id):
    

def convert_npy2cfg(population_id, n_ramdom, nqirr, Temperature):
    """
    Convert a *.npy to a cfg file.
    """
    # Load the *.npy file

    dyn = CC.Phonons.Phonons(f"dyn_pop{population_id-1}_", nqirr = nqirr)
    # Load the original ensemble (first population with 1000 configurations)
    ens = sscha.Ensemble.Ensemble(dyn, Temperature, dyn.GetSupercell())
    ens.load_bin(f"ensembles", population_id = population_id)
    ens.save("ensemble_dat", population = population_id, N = n_ramdom)
    convert_dat2cfg(population_id, n_ramdom)


#-----------------------0.frequently changed parameters-----------------------#
ONSET_DYN_POP_IDX = 15
MAX_POPULATION = 30
TARGET_PRESSURE = 200 # GPa
NQIRR = 4

#-----------------------0.frequently changed parameters-----------------------#




for POPULATION in range(ONSET_DYN_POP_IDX, MAX_POPULATION+1):
    
    MATDYN = f"dyn_pop{str(POPULATION-1)}_"
    
    dyn = CC.Phonons.Phonons(MATDYN, NQIRR) #LOAD THE DYN MAT IN MEMORY
    dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
    dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL) #LOADS THE DYN IN THE SSCHA PROGRAM (class ens)
    ens.generate(N_RANDOM)       #GENERATES THE CONFIGURATIONS BASED ON DYN
    ens.save( ENS_FOLDER + str(POPULATION), POPULATION)  #SAVES THE CONFIGURATIONS ON FILE

    convert_dat2cfg(POPULATION)
    mlp calc-efs pot.mtp train.cfg calculated.cfg
    convert_cfg2dat(POPULATION)
    
    ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION, N=N_RANDOM)   #THIS IS TO USE IF YOU ALREADY HAVE THE CONFIGURATIONS AND DON'T WANT TO GENERATE NEW ONES
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ens) #LOADS THE ROUTINE FOR THE MINIMIZATION
    
    minimizer.neglect_symmetries = False
    minimizer.min_step_dyn = 0.01 #Values around 1 are good
    minimizer.min_step_struc = 0.01
    minimizer.kong_liu_ratio = 0.5 # Usually 0.5 is a good value
    minimizer.meaningful_factor = 0.0000001
    minimizer.max_ka = MAX_KA # This is the maximum number of steps (if negative = infinity)
    
    relax = sscha.Relax.SSCHA(minimizer, ase_calculator = None,
                                N_configs = N_RANDOM,
                                max_pop = 1,
                                save_ensemble = True,
                                cluster = None)
    IO_freq = sscha.Utilities.IOInfo()
    IO_freq.SetupSaving(f"minimization_data_{POPULATION}")
    relax.setup_custom_functions(custom_function_post = IO_freq.CFP_SaveAll)

    # relax.vc_relax(target_press = TARGET_PRESS, static_bulk_modulus = 300, ensemble_loc = "ensembles", restart_from_ens = True, start_pop = DYN_POP_IDX+1) # start_pop = DYN_POP_IDX+1 将迭代出的DYN_POP_IDX+1代的动力学矩阵保存下来，编号为DYN_POP_IDX+1代
    relax.relax(get_stress=True, ensemble_loc = ENS_FOLDER, restart_from_ens = True, start_pop = POPULATION) # start_pop = DYN_POP_IDX+1 将迭代出的DYN_POP_IDX+1代的动力学矩阵保存下来，编号为DYN_POP_IDX+1代
    # relax.relax在执行过程中自动保存新迭代出的dyn文件，并保存为文件名为dyn_pop*_*
    relax.minim.finalize(verbose = 2)

    relax.minim.plot_results(save_filename = f"save_filename_{POPULATION}.dat",plot = False)

    # FUNCTION TO SAVE DATA###############################
    #--------------------------5.prepare sscha minimizer--------------------------#
    print("5.prepare sscha minimizer")
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    # We set up the minimization parameters
    minim.min_step_dyn = 0.01     # The minimization step on the dynamical matrix
    minim.min_step_struc = 0.01   # The minimization step on the structure
    minim.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
    minim.gradi_op = "all" # Check the stopping condition on both gradients
    minim.meaningful_factor = 0.00001 # How much small the gradient should be before I stop?
    minim.neglect_symmetries = False # whether neglect symmetries or not when relaxing structures

    number_configs = my_hpc.batch_size * my_hpc.job_number # Total random configurations 
    max_population = MAX_POPULATION # If sscha calculation is executed to dyn_pop20_, so the calculation will be stopped.
    print("number_configs",number_configs)
    relax = sscha.Relax.SSCHA(minim, ase_calculator = espresso_calc,
                            N_configs = number_configs,
                            max_pop = max_population,
                            save_ensemble = True,
                            cluster = my_hpc)
    #--------------------------prepare sscha minimizer--------------------------#




