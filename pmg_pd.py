#!/usr/bin/env python
# -*- coding: utf-8 -*-

# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import sys
from numbers import Number
import matplotlib
matplotlib.use("pdf")

from pymatgen.analysis.phase_diagram import *
from pymatgen.core.periodic_table import Element, DummySpecie
#from pymatgen.core.composition import Composition


if __name__ == '__main__':
    infile = sys.argv[1]
    if os.path.splitext(infile)[1].lower() == 'csv':
        (elements, entries) = PDEntry.from_csv(infile) 
    else :
        entries = []
        with open(infile) as fr :
            for line in fr :
                try :
                    comp, ene = line.split(',')
                except :
                    comp, ene = line.split()
                entries.append(PDEntry(Composition(comp), float(ene)))

    pd = PhaseDiagram(entries)
    stable_formulas = [ent.name for ent in pd.stable_entries]
    # for formula in expected_stable:
    for formula in stable_formulas:
        print (formula)
    plotter = PDPlotter(pd, show_unstable=True)
    # plotter.show()
    plotter.write_image('pd.png', image_format='png')
    for e in entries:
        ehull = pd.get_e_above_hull(e)
        print (e.composition, ehull)
