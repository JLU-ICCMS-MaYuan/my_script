from __future__ import print_function
import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.calculators

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Cluster
import spglib

from numpy import *
import numpy as np

import os

#-----------------------0.frequently changed parameters-----------------------#
BATCH_SIZE = 20     # We can also choose in how many batch of jobs we want to submit simultaneously
JOB_NUMBER = 5     # and how many configurations for each job. # In this way we submit 25 jobs, each one with 8 configurations (overall 200 configuration at time)
ONSET_DYN_POP_IDX = 0
MAX_POPULATION = 10
TARGET_PRESSURE = 200 # GPa
#-----------------------0.frequently changed parameters-----------------------#




#------------------------1.Prepare QE input parameters------------------------#
print("1.Prepare QE input parameters")
atom_masses = {
  "Ce": 140.116,
  "Sc": 44.955912,
   "H": 1.00794
         }

bool_value = True
pseudo = {"Ce": "Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf",
          "Sc": "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF",
          "H" : "H.pbe-kjpaw_psl.1.0.0.UPF"}

input_data = {
               'control' : {
                   'disk_io' : 'none',
                   'pseudo_dir' : '/work/home/mayuan/workplace/5.calypso/35.Ce-Sc-H/4.detailed-compute/1.CeSc2H24/200GPa/0.prepare-init-dyn/pp',
                   'tstress' : bool_value,
                   'tprnfor ': bool_value
                           },
               'system' : {
                          # Specify the basis set cutoffs
                            'ecutwfc' : 80, # Cutoff for wavefunction
                          # Information about smearing (it is a metal)
                            'occupations' : 'smearing',
                            'smearing' : 'mp',
                            'degauss' : 0.02,
                            # 'nbnd' : 100,
                          },
               'electrons' : {
                             'conv_thr' : 1e-8 ,
                             'mixing_beta' : 0.5
                             },
             }
k_points = (4,4,2) # The k points grid (you can alternatively specify a kspacing)
# k_points = "gamma" # It is suitable for very big supercell
k_offset = (0,0,0) # The offset of the grid (can increase convergence)

espresso_calc = CC.calculators.Espresso(input_data,
                                     pseudo,
                                     kpts = k_points,
                                     masses=atom_masses,
                                     koffset = k_offset)
#------------------------Prepare QE input parameters------------------------#




#-------------------------2.Load the dynamical matrix-------------------------#
print("2.Load the dynamical matrix")
# Begin from dyn_pop_idx
onset_dyn_pop_idx = ONSET_DYN_POP_IDX
irr_idx = 4
dyn = CC.Phonons.Phonons(f"dyn_pop{onset_dyn_pop_idx}_", nqirr = irr_idx)

# Apply the sum rule and delete the imaginary modes
dyn.Symmetrize()
dyn.ForcePositiveDefinite()
#-------------------------Load the dynamical matrix-------------------------#




#-----------------------3.Prepare random configurations-----------------------#
print("3.Prepare random configurations")
# We prepare the ensemble for random configurations, namely, generate random configurations.
ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell = dyn.GetSupercell())
#-----------------------Prepare random configurations-----------------------#




#---------------------------4.Rrepare the cluster----------------------------#
print("4.Rrepare the cluster")
my_hpc = sscha.Cluster.Cluster(mpi_cmd=r"mpirun -n 48", AlreadyInCluster=True) # AlreadyInCluster=True代表本地提交任务, 既可以把当前python脚本提交到节点上, 也可以在主节点运行当前脚本 AlreadyInCluster=False代表通过ssh连接远程机器提交任务
my_hpc.workdir = os.path.join(os.getcwd(), "run_calculation")
my_hpc.binary  = "/work/home/mayuan/software/qe-7.1/bin/pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
# Then we need to specify if some modules must be loaded in the submission script. # All these information are independent from the calculation. Now we need some more specific info, like the number of processors, pools and other stuff
my_hpc.load_modules = """#!/bin/sh
#SBATCH  --job-name=qe
#SBATCH  --output=slurm-%j.out
#SBATCH  --error=slurm-%j.err
#SBATCH  --partition=liuhanyu
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1                         
#SBATCH  --exclude=node98

source /work/home/mayuan/intel/oneapi/setvars.sh
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0

# export MKL_DEBUG_CPU_TYPE=5
# export I_MPI_ADJUST_REDUCE=3

"""
my_hpc.n_cpu = 48  # We will use 32 processors
my_hpc.n_nodes = 1 # In 1 node
my_hpc.n_pool = 8  # This is an espresso specific tool, the parallel CPU are divided in 4 pools
my_hpc.batch_size = BATCH_SIZE     # We can also choose in how many batch of jobs we want to submit simultaneously
my_hpc.job_number = JOB_NUMBER     # and how many configurations for each job. # In this way we submit 25 jobs, each one with 8 configurations (overall 200 configuration at time)
my_hpc.set_timeout(3600)   # We give 25 seconds of timeout
my_hpc.time = "48:00:00" #  We can specify the time limit for each job, 5 minutes
my_hpc.setup_workdir()
#---------------------------Rrepare the cluster ----------------------------#




#--------------------------5.prepare sscha minimizer--------------------------#
print("5.prepare sscha minimizer")
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
# We set up the minimization parameters
minim.min_step_dyn = 0.01     # The minimization step on the dynamical matrix
minim.min_step_struc = 0.01   # The minimization step on the structure
minim.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
minim.gradi_op = "all" # Check the stopping condition on both gradients
minim.meaningful_factor = 0.00001 # How much small the gradient should be before I stop?
minim.neglect_symmetries = False # whether neglect symmetries or not when relaxing structures

number_configs = my_hpc.batch_size * my_hpc.job_number # Total random configurations 
max_population = MAX_POPULATION # If sscha calculation is executed to dyn_pop20_, so the calculation will be stopped.
print("number_configs",number_configs)
relax = sscha.Relax.SSCHA(minim, ase_calculator = espresso_calc,
                         N_configs = number_configs,
                         max_pop = max_population,
                         save_ensemble = True,
                         cluster = my_hpc)
#--------------------------prepare sscha minimizer--------------------------#




#-----------------------------6.Print space group-----------------------------#
# we define a function that prints the space group during the optimization
print("6.Print space group")
space_groups = []
angles = []
paras = []
all_frequencies = []
def print_spacegroup(minim):
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

    spgroup = spglib.get_spacegroup(minim.dyn.structure.get_ase_atoms(), 0.05)
    space_groups.append(spgroup)

    # We can save them in the output at each minimization step
    f = open("space_group.dat", "w")
    f.writelines(["{}) {}\n".format(i+1, x) for i,x in enumerate(space_groups)])
    f.close()

    w, p = minim.dyn.DiagonalizeSupercell()
    all_frequencies.append(w)
    # In this way the file will be updated at each step
    np.savetxt("all_frequencies.dat", all_frequencies)

relax.setup_custom_functions(custom_function_post = print_spacegroup)
IO_freq = sscha.Utilities.IOInfo()
IO_freq.SetupSaving("minimization_data")
relax.setup_custom_functions(custom_function_post = IO_freq.CFP_SaveAll)

#-----------------------------Print space group-----------------------------#




#-----------------------------7.Run the calculation-----------------------------#
print("7.Print space group")
# In this case we fix the volume (we optimize lattice parameters)
# But you can also fixe the target pressure (as done in the commented line)
if not os.path.exists("ensembles"):
    os.mkdir("ensembles")
relax.vc_relax(target_press = TARGET_PRESSURE, static_bulk_modulus = 300, ensemble_loc = "ensembles", start_pop = onset_dyn_pop_idx+1)
relax.minim.finalize(verbose = 2)
relax.minim.plot_results(save_filename = "save_filename.dat",plot = False)
#-----------------------------Run the calculation-----------------------------#

