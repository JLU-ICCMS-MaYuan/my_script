from numpy import *
import numpy as np

import os
import sys

num = int(sys.argv[1])
directory = "run_calculation"
output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly

print("number of output",len(output_files))

energies = np.zeros(len(output_files))
id_nums = []
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


    if Flag_Ener==False :
       print("ENERGY WRONG",id_number)


sort_num = sorted(id_nums)
sort_num.append(10897654)
#print(sort_num)
j=0
for i in range(0,num):
    if i==sort_num[j]:
       j=j+1
    else:
       print("without",i)
from numpy import *
import numpy as np

import os
import sys

num = int(sys.argv[1])
directory = "run_calculation"
output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")] # We select only the output files
output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly

print("number of output",len(output_files))

energies = np.zeros(len(output_files))
id_nums = []
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


    if Flag_Ener==False :
       print("ENERGY WRONG",id_number)


sort_num = sorted(id_nums)
sort_num.append(10897654)
#print(sort_num)
j=0
for i in range(0,num):
    if i==sort_num[j]:
       j=j+1
    else:
       print("without",i)
