#!/usr/bin/env python3

info = '''Dynamical matrix file
Electron-phonon coefficients for Nb4H14
  2   18   0   7.7738457   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Basis vectors 
      1.000000000    0.000000000    0.000000000
      0.000000000    1.000000000    0.000000000
      0.000000000    0.000000000    1.000000000
           1  'Nb  '    84678.9851990652
           2  'H   '    918.681110398939
    1    1      0.7500000000      0.7500000000      0.7500000000
    2    1      0.2500000000      0.2500000000      0.7500000000
    3    1      0.2500000000      0.7500000000      0.2500000000
    4    1      0.7500000000      0.2500000000      0.2500000000
    5    2      0.6549958319      0.6549958319      0.3450041681
    6    2      0.6549958319      0.3450041681      0.6549958319
    7    2      0.3450041681      0.6549958319      0.6549958319
    8    2      0.8450041681      0.8450041681      0.1549958319
    9    2      0.1549958319      0.1549958319      0.1549958319
   10    2      0.8450041681      0.1549958319      0.8450041681
   11    2      0.1549958319      0.8450041681      0.8450041681
   12    2      0.3450041681      0.3450041681      0.3450041681
   13    2      0.5000000000      0.0000000000      0.0000000000
   14    2      0.0000000000      0.5000000000      0.0000000000
   15    2      0.0000000000      0.0000000000      0.5000000000
   16    2      0.5000000000      0.0000000000      0.5000000000
   17    2      0.0000000000      0.5000000000      0.5000000000
   18    2      0.5000000000      0.5000000000      0.0000000000

# 2 elements, 18 atoms, 7.7738457 bohr celldm(1)
# 晶格单位是alat, 对应的scffit.in scf.in中都要用alat
# 坐标单位也是alat
'''

print(info)

import os
import sys
from argparse import ArgumentParser
from collections import Counter 

def get_struct_info(inputfile: str):
    with open(inputfile, 'r') as f:
        lines = f.readlines()
    
    nelements = int(lines[2].strip().split()[0])
    
    totnatoms = int(lines[2].strip().split()[1])
    
    celldm1 = float(lines[2].strip().split()[3])  # 单位是bohr
    celldm1 = float(f"{celldm1:.10f}")
    
    cell = [list(map(float, line.strip().split())) for line in lines[4:7]]
    cell = [[float(f"{x:.10f}") for x in vector] for vector in cell]
    
    elename = {int(line.strip().split("'")[0].strip()):line.strip().split("'")[1].strip() 
               for line in lines[7:7+nelements]} # elename = {1: 'Nb', 2: 'H'}
    
    eleweights = {line.strip().split("'")[1].strip():float(f'''{line.strip().split("'")[2].strip()}''') 
                  for line in lines[7:7+nelements]}
    
    coords = [list(map(float, line.strip().split()[1:])) for line in lines[7+nelements:7+nelements+totnatoms]]
    coords = [[float(f"{c:.10f}") for c in coord] for coord in coords]

    natoms = [int(line.strip().split()[1]) for line in lines[7+nelements:7+nelements+totnatoms]]
    element_counts = Counter(natoms)  # element_counts = {2:14, 1:4}
    # return nelements, totnatoms, celldm1, cell, elename, eleweights, coords
    return nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords

def output_vasp(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords):
    with open('POSCAR', 'w') as f:
        f.write('TITEL\n')
        f.write('{:<.10f}\n'.format(celldm1*0.529177210903))
        for vector in cell:
            f.write('{:<.10f} {:<.10f} {:<.10f}\n'.format(vector[0], vector[1], vector[2]))
        for num, name in elename.items():
            f.write('{}  '.format(name))
        else:
            f.write('\n')
        for num, name in elename.items():
            f.write('{}  '.format(element_counts[num]))
        else:
            f.write('\n')
        f.write('Direct\n')
        for coord in coords:
            f.write('{:<.10f}   {:<.10f}   {:<.10f}\n'.format(coord[1], coord[2], coord[3]))

def output_qe(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords):

    scf_fit_in  = "scffit.in"
    system_name = ''.join([ name+str(element_counts[num]) for num, name in elename.items()]) # element_counts = {2:14, 1:4} elename = {1: 'Nb', 2: 'H'}
    print(system_name)
    workpath_pppath = input("please input pp path with absolute path")
    with open(scf_fit_in, "w") as qe:
        qe.write("&CONTROL\n")
        qe.write(" calculation='scf',              \n")
        qe.write(" restart_mode='from_scratch',    \n")
        qe.write(" prefix='{}',                    \n".format(system_name))
        qe.write(" pseudo_dir='{}',                \n".format(workpath_pppath)
        qe.write(" verbosity = 'high',             \n")  
        qe.write(" outdir='./tmp',                 \n")
        qe.write(" forc_conv_thr = 1.0d-6,         \n")
        qe.write(" etot_conv_thr = 1.0d-7,         \n")
        qe.write(" tstress=.true.,                 \n")
        qe.write(" tprnfor=.true.,                 \n")

        qe.write("/\n")

        qe.write("&SYSTEM\n")
        qe.write(" ibrav=0,                        \n")  # 设置ibrav=0，这时需要在输入文件中写入CELL_PARAMETERS，即CELL的基矢量. alat bohr angstrom alat 由 celldm(1)或A定义的晶格常数单位
        qe.write(" nat={},                         \n".format(totnatoms))
        qe.write(" ntyp={},                        \n".format(len(element_counts))))
        qe.write(" occupations = 'smearing',       \n")
        qe.write(" smearing = 'gauss',             \n")
        qe.write(" degauss = 0.02,                 \n")
        qe.write(" ecutwfc = 80,                   \n")
        qe.write(" ecutrho = 960,                  \n")
        qe.write(" lspinorb = .false.,             \n")
        qe.write(" noncolin = .false.,             \n")
        qe.write(" la2F = .true.,                  \n")
        qe.write(" celldm(1) = {},                 \n".format(celldm1))
        qe.write("/\n")

        qe.write("&ELECTRONS\n")
        qe.write(" conv_thr = 1.0d-9,              \n")
        qe.write(" mixing_beta = 0.7,              \n")
        qe.write(" electron_maxstep = 200,         \n")
        qe.write("/\n")

        qe.write("ATOMIC_SPECIES                   \n")
        for num, name in elename.items():
            qe.write(" {:<5}  {:<10}  {:<50} \n".format(name, "unknown", name+'.UPF'))
        qe.write(r"CELL_PARAMETERS {alat}           \n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
        for cell_p in cell:
            qe.write("{}\n".format(cell_p))
        qe.write(r"ATOMIC_POSITIONS {alat}          \n")
        for coord in coords:
            qe.write('{:<4}   {:<.10f}   {:<.10f}   {:<.10f}\n'.format(elename[int(coord[0])], coord[1], coord[2], coord[3]))
        qe.write("K_POINTS {automatic}             \n")
        qe.write(" 8 8 8 0 0 0                     \n")




if __name__ == '__main__':
    parser = ArgumentParser(description="Convert QE dynamical matrix file to POSCAR format")
    parser.add_argument('-i', '--inputfile', required=True, help="Input file name")
    parser.add_argument('-o', '--output_format', choices=['qe', 'vasp'], default='vasp', help="Output format (default: vasp)")
    
    args = parser.parse_args()
    inputfile = args.inputfile
    output_format = args.output_format

    nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords = get_struct_info(inputfile)
    if output_format == 'vasp':
        output_vasp(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords)
    elif output_format == 'qe':
        output_qe(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords)
