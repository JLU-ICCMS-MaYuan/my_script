from numpy import *
import numpy as np

import os

# Ok now we will use python to to parse the ensemble and generate input files for the ab-initio run

# First of all, we must prepare a generic header for the calculation
typical_espresso_header = """
&control
    calculation = "scf"
    tstress = .true.
    tprnfor = .true.
    disk_io = "none"
    pseudo_dir = "/lustre/home/h240012/work/SSCHA/0_phon/pp"
/
&system
        ecutwfc = 70
        occupations = "smearing"
        smearing = "mp"
        degauss = 0.02
        nat = 216
        ntyp = 3
        ibrav = 0
/
&electrons
    conv_thr = 1d-8
    mixing_beta = 0.5
/
ATOMIC_SPECIES
 Ce     140.116     Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf
 Sc     44.955912   Sc.pbe-spn-kjpaw_psl.1.0.0.UPF
 H      1.00794     H.pbe-kjpaw_psl.1.0.0.UPF
K_POINTS {automatic}
2 2 3 0 0 0
"""

# Now we need to read the scf files
all_scf_files = [os.path.join("data_ensemble_manual", f) for f in os.listdir("data_ensemble_manual") if f.startswith("scf_population12_")]


# We will generate the input file in a new directory
if not os.path.exists("run_calculation"):
    os.mkdir("run_calculation")

for file in all_scf_files:
    # Now we are cycling on the scf_ files we found.
    # We must extract the number of the file
    # The file is the string "data_ensemble_manual/scf_population1_X.dat"
    # Therefore the X number is after the last "_" and before the "." character
    # We can split before the string file at each "_", isolate the last part "X.dat"
    # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
    # then we convert the "X" string into an integer
    number = int(file.split("_")[-1].split(".")[0])

    # We decide the filename for the espresso input
    # We will call it run_calculation/espresso_run_X.pwi
    filename = os.path.join("run_calculation", "espresso_run_{}.pwi".format(number))

    # We start writing the file
    with open(filename, "w") as f:
        # We write the header
        f.write(typical_espresso_header)

        # Load the scf_population_X.dat file
        ff = open(file, "r")
        structure_lines = ff.readlines()
        ff.close()

        # Write the content on the espresso_run_X.pwi file
        # Note in the files we specify the units for both the cell and the structure [Angstrom]
        f.writelines(structure_lines)


