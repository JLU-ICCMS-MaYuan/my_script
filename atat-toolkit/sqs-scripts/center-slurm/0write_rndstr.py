#!/usr/bin/env python3
import sys

from pymatgen.core.structure import Structure

struct = Structure.from_file(sys.argv[1])
struct.to(filename='rndstr.in')