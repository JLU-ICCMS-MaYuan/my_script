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

IRR_IDX = 4
DYN_POP_IDX = 5 # 准备读入的dyn的代编号
N_CONFIGS = 100
TARGET_PRESSURE = 300

dyn = CC.Phonons.Phonons(f"dyn_pop{DYN_POP_IDX}_", nqirr = IRR_IDX)

dyn.Symmetrize()
dyn.ForcePositiveDefinite()

ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell = dyn.GetSupercell())
#ensemble.load("ensembles", population = 4, N = 100)
ensemble.load_bin("ensembles", population_id = DYN_POP_IDX+1) # 读取DYN_POP_IDX+1代的受力，能量，结构，应力
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

minim.min_step_dyn = 0.01         # The minimization step on the dynamical matrix
minim.min_step_struc = 0.01       # The minimization step on the structure
minim.kong_liu_ratio = 0.5        # The parameter that estimates whether the ensemble is still good
minim.gradi_op = "all"            # Check the stopping condition on both gradients
minim.meaningful_factor = 0.00001 # How much small the gradient should be before I stop?
minim.neglect_symmetries = False  # whether neglect symmetries or not when relaxing structures

relax = sscha.Relax.SSCHA(minim, ase_calculator = None,
                         N_configs = N_CONFIGS,
                         max_pop = 1,
                         save_ensemble = True,
                         cluster = None)
IO_freq = sscha.Utilities.IOInfo()
IO_freq.SetupSaving(f"minimization_data_{DYN_POP_IDX+1}")

relax.setup_custom_functions(custom_function_post = IO_freq.CFP_SaveAll)

relax.vc_relax(target_press = TARGET_PRESSURE, static_bulk_modulus = 300, ensemble_loc = "ensembles", restart_from_ens = True, start_pop = DYN_POP_IDX+1) # start_pop = DYN_POP_IDX+1 将迭代出的DYN_POP_IDX+1代的动力学矩阵保存下来，编号为DYN_POP_IDX+1代
relax.minim.finalize(verbose = 2)

relax.minim.plot_results(save_filename = f"save_filename_{DYN_POP_IDX+1}.dat",plot = False)
