#!/bin/bash
from ase.eos import EquationOfState

with open('eos.in', 'r') as f:
    lines = f.readlines()

v_E = lines[5:]
volumes  = [float(v.split()[0]) for v in lines[5:]]
energies = [float(e.split()[1]) for e in lines[5:]]
            
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
eos.plot('eos.png')

