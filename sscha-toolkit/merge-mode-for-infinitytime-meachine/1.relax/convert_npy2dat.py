import os

import cellconstructor as CC
from sscha.Ensemble import Ensemble

NQIRR = 4
DYN_POP_IDX = 5  # 准备读入的dyn的代编号
END_POP_IDX=DYN_POP_IDX+1

dyn = CC.Phonons.Phonons(f"dyn_pop{DYN_POP_IDX}_", nqirr = NQIRR)               
ens = Ensemble(dyn0=dyn, T0=0, supercell=dyn.GetSupercell())

ens.load_bin(data_dir="ensembles", population_id = END_POP_IDX)
ens.save("data_ensemble_manual", population = END_POP_IDX)


