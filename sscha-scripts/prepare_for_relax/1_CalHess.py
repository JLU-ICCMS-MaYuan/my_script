from __future__ import print_function
from __future__ import division

# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons

# Import the SCHA modules
import sscha, sscha.Ensemble


# Here the input information
DATA_DIR ="data_ensemble_manual"
N_RANDOM = 1000
DYN_PREFIX =  'dyn_pop11_'  # dyn mat that generated the last population
FINAL_DYN =   'dyn_pop12_'    # SSCHA dyn mat obtained with the last minimization 
SAVE_PREFIX = 'V3_Hessian.dyn'  # Free energy Hessian dynamical matrices
NQIRR = 4
Tg = 0
T =  0
POPULATION = 12
INCLUDE_V4 = False # True to include the 4th-order SSCHA FC term to calculate the Hessian 
print("Loading the original dynamical matrix...")
dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
print("Loading the current dynamical matrix...")
final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)

print("Loading the ensemble...")
# Now we load the ensemble
ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
#ens.load_bin(DATA_DIR, population_id = POPULATION)
ens.load(DATA_DIR, population = POPULATION , N=N_RANDOM )

print("Updating the importance sampling...")
ens.update_weights(final_dyn, T)

print("Computing the free energy hessian...")
# Set get_full_hessian to false to have only the odd correction
# Usefull if you want to study the convergence with the number of configuration
dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)

print("Saving the hessian to {}...".format(SAVE_PREFIX))
dyn_hessian.save_qe(SAVE_PREFIX)
print("Done.")

