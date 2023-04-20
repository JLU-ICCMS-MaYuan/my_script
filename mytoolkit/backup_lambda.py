#!/usr/bin/env python3

import os
import shutil
import sys

mu = sys.argv[1]
print("你使用的电荷屏蔽常数是：", mu)


if os.path.exists('lambda.in'):
    shutil.copy('lambda.in', mu+'-lambda.in')
    print("lambda.in rename finish")
else:
    print("lambda.in doesn't exist!")


if os.path.exists('lambda.out'):
    shutil.copy('lambda.out', mu+'-lambda.out')
    print("lambda.out rename finish")
else:
    print("lambda.out doesn't exist!")


if os.path.exists('alpha2F.dat'):
    shutil.copy('alpha2F.dat', mu+'-alpha2F.dat')
    print("alpha2F.dat rename finish")
else:
    print("alpha2F.dat doesn't exist!")


if os.path.exists('INPUT'):
    shutil.copy('INPUT', mu+'-INPUT')
    print("INPUT rename finish")
else:
    print("INPUT doesn't exist!")


if os.path.exists('ALPHA2F.OUT'):
    shutil.copy('ALPHA2F.OUT', mu+'-ALPHA2F.OUT')
    print("ALPHA2F.OUT rename finish")
else:
    print("ALPHA2F.OUT doesn't exist!")


if os.path.exists('ELIASHBERG.OUT'):
    shutil.copy('ELIASHBERG.OUT', mu+'-ELIASHBERG.OUT')
    print("ELIASHBERG.OUT rename finish")
else:
    print("ELIASHBERG.OUT doesn't exist!")


if os.path.exists('ELIASHBERG_IA.OUT'):
    shutil.copy('ELIASHBERG_IA.OUT', mu+'-ELIASHBERG_IA.OUT')
    print("ELIASHBERG_IA.OUT rename finish")
else:
    print("ELIASHBERG_IA.OUT doesn't exist!")


if os.path.exists('ELIASHBERG_GAP_T.OUT'):
    shutil.copy('ELIASHBERG_GAP_T.OUT', mu+'-ELIASHBERG_GAP_T.OUT')
    print("ELIASHBERG_GAP_T.OUT rename finish")
else:
    print("ELIASHBERG_GAP_T.OUT doesn't exist!")