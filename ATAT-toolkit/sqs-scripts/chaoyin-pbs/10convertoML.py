import numpy as np
from sys import argv
from ase.atoms import Atoms
from ase.io import read

def write2my(file_path, ene_i, lat_i, ele_i, coo_i, foc_i, vir_i=None):
    lat_i = lat_i.reshape(3, 3)
    coo_i = coo_i.reshape(-1, 3)
    foc_i = foc_i.reshape(-1, 3)
    with open(file_path, 'w') as file:
        file.write(f"# total energy = {ene_i} eV\n\n")
        if vir_i is not None:
            vir_i = np.array(vir_i).reshape(-1)
            file.write(f"VIRIAL\n")
            for i in vir_i:
                file.write(f'{i:20.8f}')
        file.write("\n")

        file.write("CRYSTAL\n")
        file.write("PRIMVEC\n")

        for j in lat_i:
            for k in j:
                file.write(f'{k:20.8f}')
            file.write('\n')

        file.write("PRIMCOORD\n")
        file.write(f"{len(coo_i)} 1\n")

        for j in range(len(coo_i)):
            file.write(f'{ele_i[j]:2}')
            # coo
            for k in coo_i[j]:
                file.write(f"{k:20.8f}")
            # force
            for k in foc_i[j]:
                file.write(f"{k:20.8f}")
            file.write("\n")
    pass

def read_from_ase_atoms(atoms: Atoms):
    try:
        ene = atoms.get_potential_energy()
    except:
        ene = 0.0
    lat = atoms.get_cell()
    pos = atoms.get_positions()
    try:
        foc = atoms.get_forces()
    except:
        foc = np.zeros_like(pos)
    try:
        sts = atoms.get_stress()
        xx, yy, zz, yz, xz, xy = - sts * atoms.get_volume()
        vir = np.array(
            [[xx, xy, xz],
             [xy, yy, yz],
             [xz, yz, zz]]).reshape(-1)
    except:
        vir = np.zeros([3, 3])
    ele = atoms.get_chemical_symbols()
    return ene, lat, pos, foc, vir, ele


if __name__ == '__main__':
    print(argv[1], argv[2])
    write_dir = argv[2]
    datas = read(argv[1], index=":")
    for idx, data in enumerate(datas):
        ene, lat, pos, foc, vir, ele = read_from_ase_atoms(data)
        single_write_dir = "step" + (str(idx+1)) + '_' + write_dir 
        write2my(
            single_write_dir,
            ene_i=ene,
            lat_i=lat,
            ele_i=ele,
            coo_i=pos,
            foc_i=foc,
            vir_i=vir)
