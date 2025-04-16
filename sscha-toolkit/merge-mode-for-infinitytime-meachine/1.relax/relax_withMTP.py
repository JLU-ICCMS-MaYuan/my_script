from __future__ import print_function

import cellconstructor as CC

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import spglib

from numpy import *
import numpy as np
from ase import Atoms
from ase.units import GPa

import os
import sys
import time
import shutil
from pathlib import Path 
from collections import defaultdict

angstrom2bohr = 1.88972612462577
eV2Ry = 0.0734986443513116
Ry2eV = 1.0 / eV2Ry
bohr2angstrom = 1.0 / angstrom2bohr

#-----------------------0.frequently changed parameters-----------------------#
# They are often modified
ONSET_DYN_POP_IDX = 41
MAX_POPULATION = 42
N_RANDOM = 10

# They are modified only once.
NQIRR = 4
TEMPERATURE = 0
TARGET_PRESSURE = 100 # GPa
SYMBOLS = ['H']
LOCALRUN_MTP = True
#-----------------------0.frequently changed parameters-----------------------#


def find_files(directory, pw_name_type=None, numberofpw=None):
    """
    查找目录中匹配特定模式的文件
    
    参数:
        directory (str): 搜索目录
        filename (str): 文件名匹配模式（支持前缀检测）
        numberofpw (int): 若文件名以 ESP_ 或 espresso_run_ 开头，生成指定数量的文件路径
    
    返回:
        list[Path]: 匹配的文件路径列表
    """
    result = []
    
    # 如果文件名以 ESP_ 开头且指定了数量
    if pw_name_type.startswith('ESP_'):
        for num in range(numberofpw):
            pwfile = Path(directory).joinpath(f'ESP_{num}.pwi')
            result.append(pwfile)
        return result
    # 如果文件名以 espresso_run_ 开头且指定了数量
    elif pw_name_type.startswith('espresso_run_'):
        for num in range(numberofpw):
            pwfile = Path(directory).joinpath(f'espresso_run_{num}.pwi')
            result.append(pwfile)
        return result
    # 如果文件名以其他模式开头
    else:
        sys.exit(1)
        print("Error: 未提供有效的文件名")
    
def dump_cfg(frames, filename, symbol_to_type, mode='w'):
    with open(filename, mode) as f:
        for atoms in frames:
            print(atoms)
            ret = ''
            ret += 'BEGIN_CFG\n'
            ret += 'Size\n{}\n'.format(len(atoms))
            try:
                cell = atoms.get_cell()[:]
                ret += 'Supercell\n{} {} {}\n{} {} {}\n{} {} {}\n'.format(*cell[0], *cell[1], *cell[2])
            except:
                pass
            cartes = atoms.positions
            has_forces = False
            has_constraints = False
            try:
                atoms.info['forces'] = atoms.get_forces()
            except:
                pass
            fields = ['id', 'type', 'cartes_x', 'cartes_y', 'cartes_z']
            if 'forces' in atoms.info:
                fields.extend(['fx', 'fy', 'fz'])
                forces = atoms.info['forces']
                has_forces = True

            flags = np.array(['T']*len(atoms))
            if atoms.constraints:
                atoms.info.append('mvable')
                for constr in atoms.constraints:
                    flags[constr.index] = 'F'
                has_constraints = True

            ret += 'AtomData: ' + ' '.join(fields) + '\n'
            for i, atom in enumerate(atoms):
                atom_info = '{} {} {} {} {} '.format(i + 1, symbol_to_type[atom.symbol], *cartes[i])
                if has_forces:
                    atom_info += '{} {} {} '.format(*forces[i])
                if has_constraints:
                    atom_info += flags[i]
                ret += atom_info + '\n'
            try:
                atoms.info['energy'] = atoms.get_potential_energy()
            except:
                pass
            if 'energy' in atoms.info:
                ret += 'Energy\n{}\n'.format(atoms.info['energy'])
            try:
                atoms.info['stress'] = atoms.get_stress()
            except:
                pass
            if 'stress' in atoms.info:
                stress = atoms.info['stress'] * atoms.get_volume() * -1.
                ret += 'PlusStress: xx yy zz yz xz xy\n{} {} {} {} {} {}\n'.format(*stress)
            if 'identification' in atoms.info:
                ret += 'Feature identification {}\n'.format(atoms.info['identification'])
            ret += 'END_CFG\n'
            f.write(ret)

def convert_dat2cfg(data_dir, symbols, population, n_random):
    
    type_to_symbol = {i: j for i, j in enumerate(symbols)}  # i: [0,1,2], j: [Ce, Sr, H]
    symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]

    cfg_file = f"init_{str(population)}.cfg"
    cfg_file = open(cfg_file, "w")
    FORCE_FILES = True
    energy_file_name = Path(data_dir).joinpath(f"energies_supercell_population{population}.dat")
    with open(energy_file_name, "r") as f:
        energy_lines = [l.strip() for l in f.readlines()]

    assert n_random == len(energy_lines), "Number of random structures does not match the number of energy lines."

    for j in range(n_random):
        structure_file_name = Path(data_dir).joinpath(f"scf_population{population}_{j+1}.dat")
        with open(structure_file_name, "r") as f:
            lines = [l.strip() for l in f.readlines()]

        try:
            force_file_name = Path(data_dir).joinpath(f"forces_population{population}_{j+1}.dat")
            with open(force_file_name, "r") as f:
                force_lines = [l.strip() for l in f.readlines()]
        except: 
            FORCE_FILES = False

        try:
            pressure_file_name = Path(data_dir).joinpath(f"pressures_population{population}_{j+1}.dat")
            with open(pressure_file_name, "r") as f:
                pressure_lines = [l.strip() for l in f.readlines()]
        except:
            pressure_lines = [
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            ]

        SIZE = len(lines)-6
        cell = np.zeros((3, 3))
        atoms = np.zeros((SIZE,3))
        forces = np.zeros((SIZE,3))
        atm_type = [None] * SIZE
        pressures = np.zeros((3, 3))

        for i in range(0,3):
            cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]

        for i in range(0,SIZE):
            atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
            atm_type[i] =  lines[i+6].split()[0]
            if FORCE_FILES == True:
                forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                forces[i,:] = forces[i, :] * Ry2eV / bohr2angstrom

        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))

        for i in range(0,3):

            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]
            pressures[i, :] = pressures[i, :] * Ry2eV / bohr2angstrom / bohr2angstrom / bohr2angstrom* Volume

        cfg_file.write("BEGIN_CFG\n")
        cfg_file.write(" Size\n")
        cfg_file.write("{: >5}".format(SIZE) +"\n")
        cfg_file.write(" Supercell\n")
        for row in cell:
            cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
        cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for i in range(0,len(atoms)):
            if FORCE_FILES == True:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(symbol_to_type[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
            else:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(symbol_to_type[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*[0,0,0]))
        cfg_file.write(" Energy\n")
        cfg_file.write("     {: >13f}".format(float(energy_lines[j])*Ry2eV) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(symbol_to_type) + "\n")
        cfg_file.write(" Feature conf_number " + str(j+1) + "\n")
        cfg_file.write(" Feature population " + str(population) + "\n")
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()

def load_cfg(population, type_to_symbol):
    frames = []
    calculated_cfg = Path("run_calculation").joinpath(f"calculated_{population}.cfg")
    with open(calculated_cfg) as f:
        line = 'chongchongchong!'
        while line:
            line = f.readline()
            if 'BEGIN_CFG' in line:
                cell = np.zeros((3, 3))

            if 'Size' in line:
                line = f.readline()
                natoms = int(line.split()[0])
                positions = np.zeros((natoms, 3))
                forces = np.zeros((natoms, 3))
                energies = np.zeros(natoms)
                symbols = ['X'] * natoms
                constraints = ['T'] * natoms

            if 'Supercell' in line: 
                for i in range(3):
                    line = f.readline()
                    for j in range(3):
                        cell[i, j] = float(line.split()[j])

            if 'AtomData' in line:
                d = defaultdict(int)
                for (i, x) in enumerate(line.split()[1:]):
                    d[x] = i

                for _ in range(natoms):
                    line = f.readline()
                    fields = line.split()
                    i = int(fields[d['id']]) - 1
                    symbols[i] = type_to_symbol[int(fields[d['type']])]
                    positions[i] = [float(fields[d[attr]]) for attr in ['cartes_x', 'cartes_y' ,'cartes_z']]
                    forces[i] = [float(fields[d[attr]]) for attr in ['fx', 'fy' ,'fz']]
                    energies[i] = float(fields[d['site_en']])

                atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)
                if d['fx'] != 0:
                    atoms.info['forces'] = forces
                if d['site_en'] != 0:
                    atoms.info['energies'] = energies

            if 'Energy' in line and 'Weight' not in line:
                line = f.readline()
                atoms.info['energy'] = float(line.split()[0])

            if 'PlusStress' in line:
                line = f.readline()
                plus_stress = np.array(list(map(float, line.split())))
                atoms.info['PlusStress'] = plus_stress
                atoms.info['stress'] = -plus_stress / atoms.get_volume()
                atoms.info['pstress'] = atoms.info['stress'] / GPa

            if 'END_CFG' in line:
                frames.append(atoms)

            if 'identification' in line:
                atoms.info['identification'] = int(line.split()[2])

    return frames

def convert_cfg2dat(data_dir, population):

    frames = load_cfg(population=population, type_to_symbol=type_to_symbol)
  
    energies = np.zeros(len(frames))
    
    for i, atoms in enumerate(frames):
        energies[i] = atoms.info['energy']*eV2Ry
        force_filename = Path(data_dir).joinpath(f"forces_population{population}_{i+1}.dat")
        np.savetxt(force_filename, atoms.info['forces']*(eV2Ry/angstrom2bohr))
        
        pressure_filename = Path(data_dir).joinpath(f"pressures_population{population}_{i+1}.dat")
        pressure_eV_in_a_volume_Angstrom3 = atoms.info['PlusStress']
        pressure_Ry_per_volume_bohr3 = [[0,0,0],[0,0,0],[0,0,0]]
        pressure_Ry_per_volume_bohr3[0][0] = pressure_eV_in_a_volume_Angstrom3[0]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[1][1] = pressure_eV_in_a_volume_Angstrom3[1]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        pressure_Ry_per_volume_bohr3[2][2] = pressure_eV_in_a_volume_Angstrom3[2]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[1][2] = pressure_eV_in_a_volume_Angstrom3[3]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        pressure_Ry_per_volume_bohr3[0][2] = pressure_eV_in_a_volume_Angstrom3[4]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)
        pressure_Ry_per_volume_bohr3[0][1] = pressure_eV_in_a_volume_Angstrom3[5]  *  eV2Ry/(atoms.get_volume()*angstrom2bohr**3)     
        np.savetxt(pressure_filename, pressure_Ry_per_volume_bohr3)

    energy_file_name = Path(data_dir).joinpath(f"energies_supercell_population{population}.dat")
    np.savetxt(energy_file_name, energies)

def calculate_efs(population, localrun_MTP=False):
    # 读取所有的能量文件
    if localrun_MTP == True:
        if not os.path.exists("run_calculation"):
            os.mkdir("run_calculation")

        print(f"submit job for calculate init_{population}.cfg")
        shutil.copy(f"init_{population}.cfg", "run_calculation")
        cwd = os.getcwd()
        os.chdir("run_calculation")
        os.system("mpirun -np 4 mlp calc-efs pot.mtp init_{population}.cfg calculated_{population}.cfg")
        os.chdir(cwd)
    else:
    # 读取所有的能量文件
        slurm_scripts = f"""#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --job-name=calcefs
#SBATCH --output=calcefs_{population}.out


source /public/home/mayuan/intel/oneapi/setvars.sh --force    
mpirun -np 56 mlp calc-efs pot.mtp init_{population}.cfg calculated_{population}.cfg

"""
        if not os.path.exists("run_calculation"):
            os.mkdir("run_calculation")
        slurm_scripts_name = Path("run_calculation").joinpath(f"cal_efs.sh")
        with open(slurm_scripts_name, "w") as f:
            f.write(slurm_scripts)
        
        
        print(f"submit job for calculate init_{population}.cfg")
        shutil.copy(f"init_{population}.cfg", "run_calculation")
        cwd = os.getcwd()
        os.chdir("run_calculation")
        jobid = os.popen(f"sbatch cal_efs.sh").read().strip().split()[-1]
        os.chdir(cwd)
        
        run_flag = True
        while run_flag:
            squeue_output = os.popen(f"squeue -j {jobid}").readlines()
            if len(squeue_output) == 1:
                run_flag = False
                print(f"{jobid} is running.")
            else:
                print(f"{jobid} was stopped.")
                time.sleep(10) 
    
def print_spacegroup(minim):
    space_groups = []
    angles = []
    paras = []
    all_frequencies = []
    
    step_id = len(minim.__fe__)

    unit_cell = minim.dyn.structure.unit_cell
    print(" ==== UNIT_CELL ====")
    print(unit_cell)
    print("//////////////////////////////")

    print (" ==== STRUCTURE [A] ==== ")

    nat = minim.dyn.structure.N_atoms
    for i in range(nat):
        print ("%5s %16.8f%16.8f%16.8f" % (minim.dyn.structure.atoms[i],
                                          minim.dyn.structure.coords[i,0],
                                          minim.dyn.structure.coords[i,1],
                                          minim.dyn.structure.coords[i,2]))
    print ("/////////////////////////////")


    para = np.sqrt(unit_cell[0,0]*unit_cell[0,0]+unit_cell[0,1]*unit_cell[0,1]+unit_cell[0,2]*unit_cell[0,2])
    paras.append(para)
    np.savetxt("para.dat", paras)


    # Compute the three angles
    angle = np.zeros(3)
    for i in range(3):
        nexti = (i+1)%3
        otheri = (i+2)%3
        angle[otheri] = np.arccos( np.dot(unit_cell[i,:], unit_cell[nexti,:]) /
             (np.sqrt(np.dot(unit_cell[i,:], unit_cell[i,:])) *
              np.sqrt(np.dot(unit_cell[nexti,:], unit_cell[nexti,:])))) * 180 / np.pi
    angles.append(angle)
    np.savetxt("angle.dat", angles)

    atoms = minim.dyn.structure.get_ase_atoms()
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)
    
    spgroup = spglib.get_spacegroup(cell, 0.05)
    space_groups.append(spgroup)

    # We can save them in the output at each minimization step
    f = open("space_group.dat", "w")
    f.writelines(["{}) {}\n".format(i+1, x) for i,x in enumerate(space_groups)])
    f.close()

    w, p = minim.dyn.DiagonalizeSupercell()
    all_frequencies.append(w)
    # In this way the file will be updated at each step
    np.savetxt("all_frequencies.dat", all_frequencies)



type_to_symbol = {i: j for i, j in enumerate(SYMBOLS)}  # i: [0,1,2], j: [Ce, Sr, H]
symbol_to_type = {v: k for k, v in type_to_symbol.items()}   # 用于输出cfg时的映射 v: [Ce, Sr, H], k: [0,1,2]

for POPULATION in range(ONSET_DYN_POP_IDX, MAX_POPULATION):
    #-------------------------1.Load the dynamical matrix-------------------------#
    print("1.Load the dynamical matrix")
    dyn_pop = f"dyn_pop{str(POPULATION)}_"
    dyn = CC.Phonons.Phonons(dyn_pop, NQIRR) #LOAD THE DYN MAT IN MEMORY
    dyn.Symmetrize()                       #IN FIRST STEP ONLY: APPLIES SUM RULE
    dyn.ForcePositiveDefinite()            #IN FIRST STEP ONLY: FLIPS NEGATIVE FREQUENCIES
    #-------------------------1.Load the dynamical matrix-------------------------#



    ensembles_dirname = f"ensembles_{POPULATION+1}"
    #-------------------------2.Prepare random configurations-------------------------#
    print("2.Prepare random configurations")
    ensemble = sscha.Ensemble.Ensemble(dyn, T0 = TEMPERATURE, supercell = dyn.GetSupercell()) #LOADS THE DYN IN THE SSCHA PROGRAM (class ensemble)
    ensemble.generate(N_RANDOM)       #GENERATES THE CONFIGURATIONS BASED ON DYN
    ensemble.save(ensembles_dirname, POPULATION+1)  #SAVES THE CONFIGURATIONS ON FILE
    #-------------------------2.Prepare random configurations-------------------------#




    #-------------------------3.convert random configurations from QE to MTP-type cfg-------------------------#
    print("3.convert random configurations from QE to MTP-type cfg")
    convert_dat2cfg(data_dir=ensembles_dirname, symbols=SYMBOLS, population=POPULATION+1, n_random=N_RANDOM)
    #-------------------------3.convert random configurations from QE to MLIP cfg-------------------------#




    #-------------------------4.evalute energy force virial by MTP-------------------------#
    print("4.evalute energy force virial by MTP")
    calculate_efs(population=POPULATION+1,  localrun_MTP=LOCALRUN_MTP)
    #-------------------------4.evalute energy force virial by MTP-------------------------#




    #--------------------------5.convert cfg to sscha-type dat--------------------------#
    print("5.convert cfg to sscha-type dat")
    convert_cfg2dat(data_dir=ensembles_dirname, population=POPULATION+1)
    #--------------------------5.convert cfg to sscha-type dat--------------------------#




    #--------------------------6.load sscha-data for minimizer--------------------------#
    print("6.load sscha-data for minimizer")
    ensemble.load(data_dir=ensembles_dirname, population=POPULATION+1, N=N_RANDOM)   #THIS IS TO USE IF YOU ALREADY HAVE THE CONFIGURATIONS AND DON'T WANT TO GENERATE NEW ONES
    #--------------------------6.load sscha-data for minimizer--------------------------#




    #--------------------------7.prepare sscha minimizer--------------------------#
    print("7.prepare sscha minimizer")
    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble) #LOADS THE ROUTINE FOR THE MINIMIZATION
    # We set up the minimization parameters
    minimizer.min_step_dyn = 0.01     # The minimization step on the dynamical matrix
    minimizer.min_step_struc = 0.01   # The minimization step on the structure
    minimizer.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
    minimizer.gradi_op = "all" # Check the stopping condition on both gradients
    minimizer.meaningful_factor = 0.00001 # How much small the gradient should be before I stop?
    minimizer.neglect_symmetries = False # whether neglect symmetries or not when relaxing structures

    print("number_configs", N_RANDOM)
    relax = sscha.Relax.SSCHA(minimizer, 
                            ase_calculator = None,
                            N_configs = N_RANDOM,
                            max_pop = 1,
                            save_ensemble = True,
                            cluster = None)
    IO_freq = sscha.Utilities.IOInfo()
    relax.setup_custom_functions(custom_function_post = IO_freq.CFP_SaveAll)
    relax.setup_custom_functions(custom_function_post = print_spacegroup)
    IO_freq.SetupSaving(f"minimization_data_{POPULATION+1}")

    #--------------------------7.prepare sscha minimizer--------------------------#




    #----------------------------8.Run the calculation----------------------------#
    print("8.Run the calculation")
    relax.vc_relax(
        target_press = TARGET_PRESSURE, 
        static_bulk_modulus = 300, 
        ensemble_loc = ensembles_dirname, 
        restart_from_ens = True, 
        start_pop = POPULATION+1) # start_pop = DYN_POP_IDX+1 将迭代出的DYN_POP_IDX+1代的动力学矩阵保存下来，编号为DYN_POP_IDX+1代
    relax.minim.finalize(verbose = 2)

    relax.minim.plot_results(save_filename = f"save_filename_{POPULATION+1}.dat",plot = False)
    #----------------------------8.Run the calculation----------------------------#
