import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

import numpy as np
inter_dyn = CC.Phonons.Phonons()
dyn_sscha = CC.Phonons.Phonons("V3_Hessian.dyn",2)
dyn_coarse= CC.Phonons.Phonons("112.dyn",2)
dyn_fine  = CC.Phonons.Phonons("224.dyn",6)
inter_dyn = dyn_sscha.Interpolate(coarse_grid=[1,1,2], fine_grid=[2,2,4], support_dyn_coarse=dyn_coarse, support_dyn_fine=dyn_fine, symmetrize=True)
#dyn_222.SwapQPoints(dyn_sscha)
#dyn_222.save_qe("i")
inter_dyn.save_qe("inter_dyn_")
