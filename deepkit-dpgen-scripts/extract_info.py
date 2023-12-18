import os,numpy,glob
import numpy as np

#outcars = glob.glob('OUTCAR_*')                                                                                                                                                                   
outcars = [f'OUTCAR_{str(i)}' for i in range(800, 1450, 50)]
free_per_atom = []
virial_per_atom = []
Tforce = []

for outcar in outcars:
        natoms = int(os.popen("grep 'NIONS' %s |awk '{print $12}'"%(outcar)).read())
        free_energy = float(os.popen("grep 'energy  without' %s |awk '{print $4}'"%(outcar)).read())
        #free_energy = float(os.popen("grep 'free  energy' %s |awk '{print $5}'"%(outcar)).read())
        virial = list(map(float, os.popen("grep -A20 '\-STRESS' %s |grep Total|awk '{print $2,$3,$4,$5,$6,$7}'"%(outcar)).read().split()))
        force = list(map(float, os.popen("grep -A%d 'TOTAL-FORCE' %s |tail -n %d|awk '{print $4,$5,$6}'"%(natoms+1,outcar,natoms)).read().split()))
        print(natoms,np.array(free_energy)/natoms)#,np.array(virial)/natoms,force)
        free_per_atom.append(np.array(free_energy)/natoms)
        virial_per_atom.append(np.array(virial)/natoms)
        Tforce.append(force)

print('free_energy','force','virial')
for idx, ff in enumerate(free_per_atom):
    detla_fe = np.abs(ff - free_per_atom[-1])
    detla_force = np.abs(np.array(Tforce[idx]) - np.array(Tforce[-1]))
    detla_vir = np.abs(np.array(virial_per_atom[idx]) - np.array(virial_per_atom[-1]) )
    print(outcars[idx].split('_')[-1], np.nanmax(detla_fe)*1000, np.nanmax(detla_force)*1000, np.nanmax(detla_vir)*1000)