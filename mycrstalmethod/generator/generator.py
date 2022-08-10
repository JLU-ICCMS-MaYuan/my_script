'''
structure generator file:
content
1. base class of structure generator
2. public function
'''
import logging

from .CrystalSpecifyWyckoffs.CrystalSpecifyWyckoffs import CrystalSpecifyWyckoffs
from .CrystalStrucGenerator.CrystalStrucGenerator import CrystalStrucGenerator


logger = logging.getLogger("GENERATOR")


def get_generator(struc_type: str):
    supported = ['crystal', 'twod', 'dummy', 'wyckoffgen']
    if struc_type not in supported:
        raise ValueError(f"Unknown structure type: {struc_type}! Should be {supported}")
    elif struc_type == 'crystal':
        Generator = CrystalStrucGenerator
    elif struc_type == 'wyckoffgen':
        Generator = CrystalSpecifyWyckoffs
    else:
        raise RuntimeError("Generator BUG! Please Report or Fix")
    logger.info("Using generator %s" % Generator.__name__)
    return Generator
