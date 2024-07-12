from __future__ import print_function
import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Cluster
import spglib

from numpy import *
import numpy as np

import os
dyn = CC.Phonons.Phonons("DYNPOP", nqirr = 4)

dyn.Symmetrize()
dyn.ForcePositiveDefinite()

ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell = dyn.GetSupercell())
ensemble.load("data_ensemble_manual", population = ENDPOP, N = ENSNUM)
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

minim.min_step_dyn = 0.01     # The minimization step on the dynamical matrix
minim.min_step_struc = 0.01   # The minimization step on the structure
minim.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
minim.gradi_op = "all" # Check the stopping condition on both gradients
minim.meaningful_factor = 0.00001 # How much small the gradient should be before I stop?
relax = sscha.Relax.SSCHA(minim, ase_calculator = None,
                         N_configs = ENSNUM,
                         max_pop = 1,
                         save_ensemble = True,
                         cluster = None)
IO_freq = sscha.Utilities.IOInfo()
IO_freq.SetupSaving("ENDPOP_freqs")

relax.setup_custom_functions(custom_function_post = IO_freq.CFP_SaveAll)

#relax.vc_relax(target_press =150, static_bulk_modulus = 300, ensemble_loc = "data_ensemble_manual",restart_from_ens = True,start_pop = ENDPOP)
relax.relax(ensemble_loc = "data_ensemble_manual",restart_from_ens = True,start_pop = ENDPOP)
relax.minim.finalize(verbose = 2)

relax.minim.plot_results(save_filename = "save_filename.dat",plot = False)
                                          
