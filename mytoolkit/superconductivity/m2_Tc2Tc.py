#!/usr/bin/env python3
from argparse import ArgumentParser

import numpy as np

parser = ArgumentParser()
parser.add_argument("--x",     dest="x",     type=float, help="待求Tc的体系的掺杂比例")
parser.add_argument("--y",     dest="y",     type=float, help="参考Tc的体系的掺杂比例")
parser.add_argument("--Tcy", dest="Tcy", type=float, help="第1个参考的Tc")
parser.add_argument("--Tc0", dest="Tc0", type=float, help="第2个参考的Tc")


args = parser.parse_args()

x = args.x
y = args.y
Tcy = args.Tcy
Tc0 = args.Tc0

Tc_x = x/y*Tcy + (1-x/y)*Tc0
print(Tc_x)