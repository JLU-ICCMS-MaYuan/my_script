#!/usr/bin/env python3

import os
import sys

sscha_out = sys.argv[1]

os.system(f'grep FC {sscha_out} > AHA_fc_results.dat')
os.system(f'grep -A3 "Ab initio average stress" {sscha_out} > AHA_ab_stress_results.dat')
os.system(f'grep -A3 "STRESS TENSOR \[GPa\]" {sscha_out} > AHA_qi_stress_results.dat')
os.system(f'grep Gibbs {sscha_out} > AHA_gibbs_results.dat')