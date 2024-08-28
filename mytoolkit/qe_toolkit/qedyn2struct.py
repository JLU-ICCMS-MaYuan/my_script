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
               for line in lines[7:7+nelements]}
    
    eleweights = {line.strip().split("'")[1].strip():float(f'''{line.strip().split("'")[2].strip()}''') 
                  for line in lines[7:7+nelements]}
    
    coords = [list(map(float, line.strip().split()[2:])) for line in lines[7+nelements:7+nelements+totnatoms]]
    coords = [[float(f"{c:.10f}") for c in coord] for coord in coords]

    natoms = [int(line.strip().split()[1]) for line in lines[7+nelements:7+nelements+totnatoms]]
    element_counts = Counter(natoms)  
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
            f.write('{:<.10f}   {:<.10f}   {:<.10f}\n'.format(coord[0], coord[1], coord[2]))

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