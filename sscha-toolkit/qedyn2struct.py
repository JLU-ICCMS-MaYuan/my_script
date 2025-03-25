#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser
from collections import Counter
from pymatgen.core.periodic_table import Element

def parse_pp_names(pp_names: list[str]) -> dict:
    """Convert --pp_names input into a dictionary."""
    pp_dict = {}
    for i in range(0, len(pp_names), 2):
        element = pp_names[i]
        pp_file = pp_names[i + 1]
        pp_dict[element] = pp_file
    return pp_dict

def get_struct_info(inputfile: str):
    with open(inputfile, 'r') as f:
        lines = f.readlines()
    
    nelements = int(lines[2].strip().split()[0])
    totnatoms = int(lines[2].strip().split()[1])
    celldm1 = float(lines[2].strip().split()[3])  # 单位是bohr
    celldm1 = float(f"{celldm1:.16f}")
    
    cell = [list(map(float, line.strip().split())) for line in lines[4:7]]
    cell = [[float(f"{x:.16f}") for x in vector] for vector in cell]
    
    elename = {int(line.strip().split("'")[0].strip()):line.strip().split("'")[1].strip() 
               for line in lines[7:7+nelements]} # elename = {1: 'Nb', 2: 'H'}
    
    eleweights = {line.strip().split("'")[1].strip():float(f'''{line.strip().split("'")[2].strip()}''') 
                  for line in lines[7:7+nelements]}
    
    coords = [list(map(float, line.strip().split()[1:])) for line in lines[7+nelements:7+nelements+totnatoms]]
    coords = [[float(f"{c:.16f}") for c in coord] for coord in coords]

    natoms = [int(line.strip().split()[1]) for line in lines[7+nelements:7+nelements+totnatoms]]
    element_counts = Counter(natoms)  # element_counts = {2:14, 1:4}
    return nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords

def output_vasp(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords):
    cell_for_angstrom = np.array(cell)* celldm1 * 0.529177210903
    coords_for_angstrom = np.array(coords)[:,-3:] * 0.529177210903 * celldm1
    with open('POSCAR', 'w') as f:
        f.write('TITEL\n')
        f.write('1.0\n')
        for vector in cell_for_angstrom:
            f.write('{:<.16f} {:<.16f} {:<.16f}\n'.format(vector[0], vector[1], vector[2]))
        for num, name in elename.items():
            f.write('{}  '.format(name))
        else:
            f.write('\n')
        for num, name in elename.items():
            f.write('{}  '.format(element_counts[num]))
        else:
            f.write('\n')
        f.write('Cartesian\n')
        for coord in coords_for_angstrom:
            f.write('{:<.16f}   {:<.16f}   {:<.16f}\n'.format(coord[0], coord[1], coord[2]))

def output_qe(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords, kd, ks, pp_path, pp_names):
    scf_fit_in  = "scffit.in"
    scf_in      = "scf.in"
    system_name = ''.join([ name+str(element_counts[num]) for num, name in elename.items()]) # element_counts = {2:14, 1:4} elename = {1: 'Nb', 2: 'H'}
    print(system_name)
    
    with open(scf_fit_in, "w") as qe:
        qe.write("&CONTROL\n")
        qe.write(" calculation='scf',              \n")
        qe.write(" restart_mode='from_scratch',    \n")
        qe.write(" prefix='{}',                    \n".format(system_name))
        qe.write(" pseudo_dir='{}',                \n".format(pp_path))
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
        qe.write(" ntyp={},                        \n".format(len(element_counts)))
        qe.write(" occupations = 'smearing',       \n")
        qe.write(" smearing = 'gauss',             \n")
        qe.write(" degauss = 0.02,                 \n")
        qe.write(" ecutwfc = 80,                   \n")
        qe.write(" ecutrho = 960,                  \n")
        qe.write(" lspinorb = .false.,             \n")
        qe.write(" noncolin = .false.,             \n")
        qe.write(" la2F = .true.,                  \n")
        qe.write(" celldm(1) = {:<.16f},           \n".format(celldm1))
        qe.write("/\n")

        qe.write("&ELECTRONS\n")
        qe.write(" conv_thr = 1.0d-9,              \n")
        qe.write(" mixing_beta = 0.7,              \n")
        qe.write(" electron_maxstep = 200,         \n")
        qe.write("/\n")

        qe.write("ATOMIC_SPECIES                   \n")
        for num, name in elename.items():
            element = Element(name)
            mass = str(element.atomic_mass).strip("amu")
            qe.write(" {:<5}  {:<10}  {:<50} \n".format(name, mass, pp_names[name]))
        qe.write(r"CELL_PARAMETERS {alat}" + "\n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
        for vector in cell:
            qe.write('{:<.16f} {:<.16f} {:<.16f}\n'.format(vector[0], vector[1], vector[2]))
        qe.write(r"ATOMIC_POSITIONS {alat}" + "\n")
        for coord in coords:
            qe.write('{:<4}   {:<.16f}   {:<.16f}   {:<.16f}\n'.format(elename[int(coord[0])], coord[1], coord[2], coord[3]))
        qe.write("K_POINTS {automatic}             \n")
        qe.write(f" {kd[0]} {kd[1]} {kd[2]} 0 0 0 \n")

    with open(scf_in, "w") as qe:
        qe.write("&CONTROL\n")
        qe.write(" calculation='scf',              \n")
        qe.write(" restart_mode='from_scratch',    \n")
        qe.write(" prefix='{}',                    \n".format( system_name))
        qe.write(" pseudo_dir='{}',                \n".format(pp_path))
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
        qe.write(" ntyp={},                        \n".format(len(element_counts)))
        qe.write(" occupations = 'smearing',       \n")
        qe.write(" smearing = 'gauss',             \n")
        qe.write(" degauss = 0.02,                 \n")
        qe.write(" ecutwfc = 80,                   \n")
        qe.write(" ecutrho = 960,                  \n")
        qe.write(" lspinorb = .false.,             \n")
        qe.write(" noncolin = .false.,             \n")
        qe.write(" celldm(1) = {:<.16f},           \n".format(celldm1))
        qe.write("/\n")

        qe.write("&ELECTRONS\n")
        qe.write(" conv_thr = 1.0d-9,              \n")
        qe.write(" mixing_beta = 0.7,              \n")
        qe.write(" electron_maxstep = 200,         \n")
        qe.write("/\n")

        qe.write("ATOMIC_SPECIES                   \n")
        for num, name in elename.items():
            element = Element(name)
            mass = str(element.atomic_mass).strip("amu")
            qe.write(" {:<5}  {:<10}  {:<50} \n".format(name, mass, pp_names[name]))
        qe.write(r"CELL_PARAMETERS {alat}" + "\n")  # 如果选择angstrom单未，原子坐标选择分数坐标，即，ATOMIC_POSITIONS (crystal), 且不设置celldm(1). 这时alat和celldm(1)设置成v1的长度
        for vector in cell:
            qe.write('{:<.16f} {:<.16f} {:<.16f}\n'.format(vector[0], vector[1], vector[2]))
        qe.write(r"ATOMIC_POSITIONS {alat}" + "\n")
        for coord in coords:
            qe.write('{:<4}   {:<.16f}   {:<.16f}   {:<.16f}\n'.format(elename[int(coord[0])], coord[1], coord[2], coord[3]))
        qe.write("K_POINTS {automatic}             \n")
        qe.write(f" {ks[0]} {ks[1]} {ks[2]} 0 0 0 \n")


if __name__ == '__main__':
    parser = ArgumentParser(description="Convert QE dynamical matrix file to POSCAR format.\n You can use this scrips by:   qedyn2struct.py -i ../V3_Hessian.dyn1 -o qe -kd 8 8 8 -ks 4 4 4 -pp /work/home/mayuan/workplace/5.calypso/35.Ce-Sc-H/4.detailed-compute/1.CeSc2H24/200GPa/0.prepare-init-dyn/pp -pn Ce Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf Sc Sc.pbe-spn-kjpaw_psl.1.0.0.UPF H H.pbe-kjpaw_psl.1.0.0.UPF")
    parser.add_argument('-i', '--inputfile', required=True, help="Input file name")
    parser.add_argument('-o', '--output_format', choices=['qe', 'vasp'], default='vasp', help="Output format (default: vasp)")
    parser.add_argument('-kd', '--kpoints_dense', nargs=3, type=int, default=[8, 8, 8], help="Dense k-point grid for scffit.in (default: 8 8 8)")
    parser.add_argument('-ks', '--kpoints_sparse', nargs=3, type=int, default=[4, 4, 4], help="Sparse k-point grid for scf.in (default: 4 4 4)")
    parser.add_argument('-pp', '--pp_path', help="Path to pseudopotential files")
    parser.add_argument('-pn', '--pp_names', nargs='+', help="Pseudopotential file names for each element")
    
    args = parser.parse_args()
    inputfile = args.inputfile
    output_format = args.output_format
    kd = args.kpoints_dense
    ks = args.kpoints_sparse
    pp_path = args.pp_path

    nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords = get_struct_info(inputfile)
    if output_format == 'vasp':
        output_vasp(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords)
    elif output_format == 'qe':
        pp_names = parse_pp_names(args.pp_names)  # 将 --pp_names 转换为字典
        output_qe(nelements, totnatoms, celldm1, cell, elename, eleweights, natoms, element_counts, coords, kd, ks, pp_path, pp_names)
        print("You have to confirm the 4 iterms, namely \n`pp path`, \n`ATOMIC_SPECIES`, \n`K_POINTS in scffit.in`, \n`K_POINTS in scf.in`.")
#        print("You can use this scrips by : qedyn2struct.py -i ../V3_Hessian.dyn1 -o qe -kd 8 8 8 -ks 4 4 4 -pp /work/home/mayuan/workplace/5.calypso/35.Ce-Sc-H/4.detailed-compute/1.CeSc2H24/200GPa/0.prepare-init-dyn/pp -pn Ce Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf Sc Sc.pbe-spn-kjpaw_psl.1.0.0.UPF H H.pbe-kjpaw_psl.1.0.0.UPF")
