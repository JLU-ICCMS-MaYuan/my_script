from pymatgen.core.structure import Structure

from structuregenerator.clathrate import HydrideClathrate

def checkclathrate(pmg_struct: Structure):
    clathrate = HydrideClathrate(pmg_struct)
    # if clathrate.remain_H_ratio > 0.44:
    #     return False
    # if clathrate.fraction_of_hydrogen_volume < 0.3:
    #     return False
    # if clathrate.shr_num_avg < 1.6:
    #     return False
    # if clathrate.cage_regularity_avg > 0.05:
    #     return False
    # if clathrate.h2h_network_regularity_avg > 0.12:
    #     return False
    if clathrate.h2h_1nearest[0][0] < 0.9:
        return False
    if clathrate.nonh2nonh_1nearest[0][0] < 2.2:
        return False
    
    return True
    