#!/usr/bin/env python3
import sys

from pymatgen.core.structure import Structure

print("If you use the `write_rndstr.py`, there have to be a *cif format file that include disorder structure information.")
print("You just need to input `write_rndstr.py FILENAME`, for example, `write_rndstr.py LaH10-disorder.cif`")
struct = Structure.from_file(sys.argv[1])
struct.to(filename='rndstr.in')

print("Write rndstr.in done !!!")
