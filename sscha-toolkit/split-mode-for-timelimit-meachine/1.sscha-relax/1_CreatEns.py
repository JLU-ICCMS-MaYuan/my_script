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
# We now load the dynamical matrix
# Let us load the starting dynamical matrix
dyn = CC.Phonons.Phonons("DYNPOP", nqirr = 4)

# Apply the sum rule and delete the imaginary modes
dyn.Symmetrize()
dyn.ForcePositiveDefinite()
# We setup an ensemble for the SSCHA at T = 100 K using the density matrix  from the dyn dynamical matrix
ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell= dyn.GetSupercell())

# We generate 10 randomly displaced structures in the supercell
ensemble.generate(N = ENSNUM)

ensemble.save("data_ensemble_manual", population = POPID)
