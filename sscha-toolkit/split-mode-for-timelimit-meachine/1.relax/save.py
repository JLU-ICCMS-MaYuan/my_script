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


dyn = CC.Phonons.Phonons("dyn_pop11_", nqirr = 4)
# Load the original ensemble (first population with 1000 configurations)
ens = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())
ens.load("data_ensemble_manual", population = 12, N = 1000)
ens.save_bin("pop12_1000",population_id = 12)

