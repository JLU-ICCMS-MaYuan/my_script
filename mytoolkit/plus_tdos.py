#!/usr/bin/env python3

import numpy as np
import pandas as pd

omega_phtdos_phpdos = pd.read_csv("phdos_proj2eles.csv")
tdos = omega_phtdos_phpdos['tdos'].values
tdos = omega_phtdos_phpdos['La'].values
tdos = omega_phtdos_phpdos['Ce'].values
tdos = omega_phtdos_phpdos['Be'].values
tdos = omega_phtdos_phpdos['H'].values
omega = omega_phtdos_phpdos['omega'].values
print(tdos)
print(omega)

domega = 4.082166131000001
print("tdos normalized")
print(np.sum(tdos*domega))

print("tdos normalized")
print(np.sum(tdos*domega))