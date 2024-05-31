#!/usr/bin/env python3
from argparse import ArgumentParser

import numpy as np

parser = ArgumentParser()
parser.add_argument("--lambda-orig", dest="lambda_orig" ,  type=float, help="已知构型的电声耦合常数")
parser.add_argument("--Nef-orig",dest="Nef_orig" ,  type=float, help="已知构型的费米能级出电子态密度")
parser.add_argument("--Tc-orig",dest="Tc_orig" ,  type=float, help="已知构型的Tc")
parser.add_argument("--Nef-x", dest="Nef_x" , type=float, help="掺杂或合金构型的费米能级出电子态密度")
parser.add_argument("--mu",dest="mu" ,  type=float, help="库伦屏蔽常数")


args = parser.parse_args()

lambda_orig = args.lambda_orig
Nef_orig= args.Nef_orig
Nef_x   = args.Nef_x
Tc_orig = args.Tc_orig
mu      = args.mu

lambda_x = lambda_orig * Nef_x / Nef_orig
print("lambda_x={}".format(lambda_x))

Tc_x = Tc_orig * np.exp(-1/(lambda_x/(1+lambda_x)-mu)) / np.exp(-1/(lambda_orig/(1+lambda_orig)-mu))

print("Tc_x={}".format(Tc_x))

# Tc_x = 169
# Tc_orig = 219
# Tc_y = 0.0625/0.125 * Tc_x + (1-0.0625/0.125)*Tc_orig
# print(Tc_y)