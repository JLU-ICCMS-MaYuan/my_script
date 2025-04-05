#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys,os
from numpy import *
import numpy as np
import time
import os
import contextlib

import ase

from ase.calculators.espresso import Espresso
from ase.visualize import view


import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Cluster

from utils_sscha_mlip.MTP_environment import MTP_Environment
from utils_sscha_mlip.Cluster_Management import Cluster_Management
import utils


#IMPORTANT PARAMETERS FOR THE WORKFLOW


GENERATE_ENS = True                                #if True, it will generate new displaced configurations. if False it will expect to find configurations already generated
T            = 0                                   #TEMPERATURE OF THE SIMULATION
N_RANDOM     = 20                                  #Number of population for each configuration
POPULATION   = 1                                   #Initial population of the simulation
POP_MAX      = 10                                  #Final population of the simulation
PREFIX       = 'PdCuH_'                            #PREFIX OF THE DYNAMICAL MATRICES
MATDYN       = PREFIX + str(POPULATION-1) + '.dyn' #NAME OF THE MATRICES
SUPERCELL    = (2,2,2)                             #size of the supercell 
NQIRR        = 6                                   # Number of irreducible matrices
ENS_FOLDER   = 'TMP_gen'                           #FOLDER FOR THE FORCES,PRESSURES,AND CONFIGURATIONS
KL_VAR       = 0.4                                 #Kong-Liu parameter (see SSCHA)
MAX_KA       = 500                                 #Maximum number of steps for the SSCHA

RELAX_FLAG          = False                        #Relaxes the cell if True
static_bulk_modulus = 10                           #Guess for the bulk modulus in GPa
target_press        = 40                           #Target pressure in GPa
fix_cell_shape      = False                        #perform an isotropic variation of the cell if True



folder     = os.getcwd() + "/SSCHA"                #Folder where to execute the dft calculations


header = utils.SCF_FILE                           #Header of the Quantum Espresso scf inputs

CALCULATOR = "ESPRESSO"                           # either "ESPRESSO" or "VASP" This flag chooses the program for the self consistent calculations
ADDRESS    = "Your cluster address "              #Your cluster address if you want to run your DFT calculations on a cluster from a local machine
marker     = "PdCuH"    #str(os.urandom(2).hex()) #Header for the scf jobs   

Submitter  = Cluster_Management(PREFIX,POPULATION,ENS_FOLDER,ADDRESS,folder, marker,N_RANDOM)   #Starting an environment for the Cluster management. Check the library for further details. By default this routine will run a local calculation. 
Submitter.base_script_VASP = utils.base_script_VASP           #Assings some slurm scripts for the scf calculations. Check utils.py file
Submitter.base_script_QE = utils.base_script_QE               #Assings some slurm scripts for the scf calculations. Check utils.py file
utils.create_scripts(folder,CALCULATOR)                       #generates folders 

TRAINED_FLAG = False                                                                   #Place true only if the SSCHA run interrupted after the MTP potential has already been trained
GAMMA        = 200                                                                     #Value of the Gamma selet for the active learning protocol
MLIP_PATH    = '/work/home/mayuan/code/mlip-2/bin/mlp'             #path to the MLIP2 execitable
MTP_handler  = MTP_Environment(POPULATION,ENS_FOLDER,GAMMA,ADDRESS,folder,MLIP_PATH)   #Prepares an handler to manage the MLIP interface

#IMPORTANT VARIABLES THAT NEED TO BE STORED THROUGHOUT THE EXECUTION OF THE SCRIPT
# energies (array of energy) contains all the dft/Mlp energies that need to be combined


for POPULATION in range(POPULATION, POP_MAX+1):

    MATDYN = f"dyn_pop{str(POPULATION-1)}_"
    Submitter.POPULATION = POPULATION
    MTP_handler.POPULATION = POPULATION
    
    #################################################################################
    dyn = CC.Phonons.Phonons(MATDYN, NQIRR) #LOAD THE DYN MAT IN MEMORY
    print(dyn.GetSupercell())               #CHECK THE SIZE OF SUPERCELL
    if POPULATION == 1:
        dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
        dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
    dyn.save_qe(f"dyn_pop{str(POPULATION-1)}_")                 #iN FIRST STEP ONLY: SAVES NEW DYN 
    
    w_s, pols = dyn.DiagonalizeSupercell()  #GETS EIGEN- VALIES AND VECTORS
    print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))
    ##################################################################################
    
    
    ##################################################################################
    dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
    dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL) #LOADS THE DYN IN THE SSCHA PROGRAM (class ens)
    if GENERATE_ENS == True:
        ens.generate(N_RANDOM)       #GENERATES THE CONFIGURATIONS BASED ON DYN
        #view(ens.structures[0].get_ase_atoms()) #VIEWS FIRST CONFIGURATION
        # Save the ensemble
        ens.save( ENS_FOLDER + str(POPULATION), POPULATION)  #SAVES THE CONFIGURATIONS ON FILE
    else:
        ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION,N=N_RANDOM,load_noncomputed_ensemble=True)
    ##################################################################################



    ##################################################################################
    #PREPARES CFG EMPTY FILES
    if MTP_handler.GAMMA > 0: MTP_handler.Generate_CFG()
    ##########################################################################################################################
    
    ##################################################################################
    #LOADS AND SUBMIT JOBS ON THE CLUSTER
    
    #CHECKS MINDIST JUST TO BE SURE NOTHING IS GOING WRONG
    if MTP_handler.GAMMA > 0: MTP_handler.Calc_GRADE()
    print("DONE")
   #####################################################################################
  
   
    
    print("creating the gamma table") 
    ###################################################################################
    if MTP_handler.GAMMA > 0: MTP_handler.Fill_GAMMA_Table()
    else: MTP_handler.Fill_GAMMA_Table_0(N_RANDOM)
    #########################################################################################
   
    
    
    ######################################################################################
    if not os.path.exists("scfin"):
        os.mkdir("scfin")
    #########################################################################################
    
    
    #########################################################################################
    #THIS CELL GENERATES THE INPUTS FOR THE QUANTUM ESPRESSO CALCULATIONS
    
    typical_espresso_header = header.format(ens.structures[0].N_atoms) 
    
    # Now we need to read the scf files
    if CALCULATOR == "ESPRESSO": Submitter.Generate_SCF(typical_espresso_header)
    elif CALCULATOR == "VASP": saved_ordering = Submitter.Generate_POSCAR()
    ##############################################################################################

    
    ##############################################################################################
    #LOADS AND SUBMIT JOBS ON THE CLUSTER
    if GENERATE_ENS == True:
        if CALCULATOR == "ESPRESSO": Submitter.Send_to_Folders(MTP_handler.conf_table_gamma)
        if CALCULATOR == "VASP": Submitter.Send_to_Folders_VASP(MTP_handler.conf_table_gamma)
        print("DONE")
    #################################################################################################
    
    print("going into queue and resubmit")
    #################################################################################################
    if CALCULATOR == "ESPRESSO": Submitter.Queue_and_resubmit(MTP_handler.conf_table_gamma)
    if CALCULATOR == "VASP": Submitter.Queue_and_resubmit_VASP(MTP_handler.conf_table_gamma)
    #######################################################################################################################################################
    
    
    #########################################################################################################################################################
    #LOADS THE forces, pressures AND energies IN THE FILES AND MEMORIES FROM THE QE INPUTS
    
    #As written, we must convert the total energy of the supercell in Ry, the forces in Ry/Bohr, and the stress in Ry/Bohr^3. Luckily quantum espresso already gives these quantities in the correct units, but be careful when using different calculators. This problem does not arise when using automatic calculators, as the SSCHA and ASE will cooperate to convert the units to the correct one. Now we will parse the Quantum ESPRESSO output looking for the energy, the forces, and the stress tensor.
    if CALCULATOR == "ESPRESSO": energies = Submitter.read_SCF()
    if CALCULATOR == "VASP": energies = Submitter.read_OUTCAR(saved_ordering)  
    print("DONE")
    ########################################################################################################################################################
    
    
    
    #########################################################################################################################################################
    #ROUTINE TO GO FROM SSCHA TMP FILE TO MLIP CFG FILE
    
    if TRAINED_FLAG == False and MTP_handler.GAMMA > 0:   
        MTP_handler.Compile_Trainingset()
        MTP_handler.Train_submit()

        print("DONE")
    
    if MTP_handler.GAMMA > 0: TOTAL_NUMBER_OF_MLIP_CONF = MTP_handler.read_CFG(energies)
    else: TOTAL_NUMBER_OF_MLIP_CONF = 0
    ######################################################################################################################################################################
    
        
   

   #######################################################################################################################################################################
    #HERE WE RUN THE MINIMIZATION
    
    
    #ens.update_weights(dyn, T=0) # CHECK THIS COMMAND ON DOCUMENTATION
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
    ens.load(data_dir=ENS_FOLDER + str(POPULATION), population=POPULATION, N=N_RANDOM)   #THIS IS TO USE IF YOU ALREADY HAVE THE CONFIGURATIONS AND DON'T WANT TO GENERATE NEW ONES
    
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ens) #LOADS THE ROUTINE FOR THE MINIMIZATION
    
    
    minimizer.neglect_symmetries = False
    minimizer.min_step_dyn = 0.01 #Values around 1 are good
    minimizer.min_step_struc = 0.01
    minimizer.kong_liu_ratio = 0.5 # Usually 0.5 is a good value
    minimizer.meaningful_factor = 0.0000001
    minimizer.max_ka = MAX_KA # This is the maximum number of steps (if negative = infinity)
    
    # FUNCTION TO SAVE DATA###############################
    all_frequencies = []
    def add_current_frequencies(minimizer):
        w, p = minimizer.dyn.DiagonalizeSupercell()
        all_frequencies.append(w)
        # In this way the file will be updated at each step
        np.savetxt("all_frequencies.dat", all_frequencies)
    #######################################################
    
    
    
    #original_stdout = sys.stdout #CHANGE STANDARD OUTPUT
    
    #with open('minimi.out', 'w') as f:
    #sys.stdout = f  #BECAUSE I LIKE TO PRINT ON FILE
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



    # HIGHLIGHTS IMPORTANT RESULTS
    print("The total free energy per unit cell is:", minimizer.get_free_energy(), " Ry")
    print("The total stress tensor is [Ry/bohr^3]:")
    print(minimizer.get_stress_tensor()[0])
    print("And the stochastic error on the stress tensor is:")
    print(minimizer.get_stress_tensor()[1])
    print("The stocastic error of the free energy instead, was:", minimizer.get_free_energy(return_error = True)[1], " Ry")
    
    #MAKES A FILE FOR THE FREQUENCIES READABLE BY XMGRACE
    res = os.system('cat "minim.freqs" | nl > "freq.dat" ')
    
    ########################################################################################################################################################################
    
    #view(minimizer.dyn.structure.get_ase_atoms()) #SHOWS FINAL STRUCTURE
    
    ###############################################################################################################################################################################
    #DOES CLEAN UP AND READIES EVERYTHING FOR THE NEXT ITERATION
    if not os.path.exists("LOG_FOLDER"):
        os.mkdir("LOG_FOLDER")
        
    if not os.path.exists("LOG_FOLDER/STEP"+str(POPULATION)):
        os.mkdir("LOG_FOLDER/STEP"+str(POPULATION))
        
    res = os.system('mv scfin minim* MLIP* dynamic.dat freq.dat TMP_gen'+ str(POPULATION) + f" dyn_pop{str(POPULATION-1)}_*  LOG_FOLDER/STEP"+str(POPULATION))