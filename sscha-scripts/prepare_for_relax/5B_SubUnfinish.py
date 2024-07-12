from numpy import *
import numpy as np
import time

import os
directory = "run_calculation"
output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly

print("number of output",len(output_files))

energies = np.zeros(len(output_files))
id_nums = []
unfinish = []
for file in output_files:
    # Get the number of the configuration.
    id_number = int(file.split("_")[-1].split(".")[0])
    id_nums.append(id_number)

    # Load the file
    ff = open(file, "r")
    lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
    ff.close()

    Flag_Ener=False
    for l in lines:
        if len(l) > 0 :
           if l.split()[0] == "!":
              Flag_Ener=True
              Flag_stress=False
              for l in lines:
                  if len(l) > 0 :
                     if l.split()[0] == "total" and l.split()[1] == "stress":
                        Flag_stress=True

              if Flag_stress==False :
                 print("Stress WRONG",id_number)
                 unfinish.append(id_number)


    if Flag_Ener==False :
       unfinish.append(id_number)
       print("ENERGY WRONG",id_number)


sort_num = sorted(id_nums)
sort_num.append(10897654)
j=0
for i in range(1,101):
    if i==sort_num[j]:
       j=j+1
    else:
       print("without",i)
       unfinish.append(i)

print(unfinish)
n_sub = len(unfinish)

#all_scf_files = [os.path.join("run_calculation", f) for f in os.listdir("run_calculation") if f.startswith("espresso_run_")]
header="""#!/bin/sh 
#------ slurm option --------#
#SBATCH --partition=mpc
#SBATCH --account=hp240139
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00

ulimit -s unlimited

#------- Program execution -------#
module load intelmpi/impi_23.2.0
module load intel/23.02.1
"""

i=50
for id_unfinish in unfinish:
    i=i+1
    filename= "un_sub_{}.sh".format(i)
    run_line="srun /home/h240012/soft/qe-7.0/bin/pw.x -nk 4 -in  espresso_run_{}.pwi > espresso_run_{}.pwo 2>&1".format(id_unfinish,id_unfinish)

    with open(filename, "w") as f:
        f.write(header)
        print(run_line,file=f)

    os.system(rf"mv {filename} run_calculation && cd run_calculation && sbatch {filename} && cd ../")
    time.sleep(3)
    if i > 100:
       break

