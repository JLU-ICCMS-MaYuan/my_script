# Script developed by H. Lambert and S. Ponce [2016]
# Updated to Python 3

import sys
import re
import numpy as np

def parse_args(args):
    extra = []
    vars = []
    current_var = None
    for arg in args:
        if arg.startswith('--'):
            current_var = arg[2:]
        else:
            if current_var is not None:
                vars.append((current_var, arg))
                current_var = None
            else:
                extra.append(arg)
    return (extra, vars)

def split_vars(vars):
    vars_values = []
    for var_name, values in vars:
        values = values.split(",")
        try:
            if any('.' in value for value in values):
                values = list(map(float, values))
            else:
                values = list(map(int, values))
        except ValueError:
            pass
        vars_values.append((var_name, values))
    vars_dict = dict(vars_values)
    return vars_dict

class FermiSurface(object):
    def __init__(self):
        self.nx = 60
        self.ny = 60
        self.nz = 60
        self.dimvec = np.array([float(self.nx), float(self.ny), float(self.nz)])
        self.fermixyz = {}
        self.gap = {}
        self.prefix = 'MgB2'
        self.nbndmin = 2
        self.nbndmax = 4

    def __repr__(self):
        return 'Fermi Surface/Mu Tensor Object'

    def cryst_to_cart(self, kvec):
        at1 = np.array([1.000000, 0.000000, 0.000000])
        at2 = np.array([-0.500000, 0.866025, 0.000000])
        at3 = np.array([0.000000, 0.000000, 1.142069])
        at = np.array([at1, at2, at3])
        outvec = np.dot(at, kvec)
        return outvec

    def pull_fermi(self, f):
        fermi_regex = re.compile(r'k\s=\s?(\-?[0-9\.]+)\s?(\-?[0-9\.]+)\s?(\-?[0-9\.]+).*?:\n\n\s+([0-9\.\-\s]+)')
        print(len(fermi_regex.findall(f)))
        for a, b, c, d in fermi_regex.findall(f):
            a = float(a)
            b = float(b)
            c = float(c)
            kvec = np.array([a, b, c])
            d = list(map(float, d.split()))
            kvec = self.cryst_to_cart(kvec)
            for i in range(len(kvec)):
                if kvec[i] < -0.001:
                    kvec[i] += 1.0
            index = [round(a) for a in np.multiply(kvec, self.dimvec)]
            d = [x - 7.4272 for x in d]
            for i in range(len(index)):
                if index[i] == 61.0:
                    index[i] = 0.0
            print(index)
            self.fermixyz[tuple(index)] = d

    def print_xsf(self, surf, title='band', band1=2, band2=3, band3=4):
        for ibnd in range(band1, band3 + 1):
            with open(f'{self.prefix}.band{ibnd}.xsf', 'w') as f1:
                print("BEGIN_BLOCK_DATAGRID_3D", file=f1)
                print(f"{self.prefix}_band_{ibnd}", file=f1)
                print(f" BEGIN_DATAGRID_3D_{self.prefix}", file=f1)
                print(f" {self.nx}  {self.ny}  {self.nz} ", file=f1)
                print("0.000000  0.000000  0.000000", file=f1)
                print("1.000000  0.577350  0.000000", file=f1)
                print("0.000000  1.154701  0.000000", file=f1)
                print("0.000000  0.000000  0.875604", file=f1)
                print("", file=f1)

                total = 0
                for z in range(self.nz):
                    for y in range(self.ny):
                        for x in range(self.nx):
                            try:
                                print(surf[x, y, z][ibnd], end="  ", file=f1)
                                total += 1
                            except TypeError:
                                print(surf[x, y, z], end="  ", file=f1)
                                print("Missing key")
                                print("0.0", end="  ", file=f1)
                            except KeyError:
                                print("Missing key")
                                print("0.0", end="  ", file=f1)
                        print("", file=f1)
                    print("", file=f1)
                print("END_DATAGRID_3D", file=f1)
                print("END_BLOCK_DATAGRID_3D", file=f1)
                print('Total number of data ', total)

if __name__ == "__main__":
    extra, vars = parse_args(sys.argv[1:])
    vars = split_vars(vars)
    print(vars, extra)

    with open(extra[0], 'r') as f:
        file_content = f.read()

    fs = FermiSurface()

    if 'fs' in vars:
        fs.pull_fermi(file_content)
        fs.print_xsf(fs.fermixyz, band1=2, band2=3, band3=4)
