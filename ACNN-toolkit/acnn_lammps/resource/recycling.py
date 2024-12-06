# acnn al recycling
#
# USAGE:
#       python ./recycling.py "root_dir" "write_dir" "fp_calculator"
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
import os
from sys import argv
from ase.io import read
import numpy as np
from ase.atoms import Atoms


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


def find_files(directory, filename):
    result = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == filename:
                result.append(os.path.join(root, file))
    return result


def get_header(path_: str):
    while True:
        if path_[0] == "." or path_[0] == "/":
            path_ = path_[1:]
        else:
            break

    return "_".join(path_.split("/"))


def get_line_index(key, lines: list, once=True):
    index_list = []
    for i, j in enumerate(lines):
        if key in j:
            index_list.append(i)
            if once:
                return index_list
    return index_list


def direct2pos(lat, pos, direction=True):
    """
    lat = [[Ax, Ay, Az],
           [Bx, By, Bz],
           [Cx, Cy, Cz]]
    Pos = [n, 3]
    :return: if True:  return cart   = [n, 3]
             if False: return direct = [n, 3]
    """
    if direction:
        return pos @ lat
    else:
        return pos @ np.linalg.inv(lat)


def read_ares_scf(path_: str):
    har2ev = 27.2113863
    bohr2a = 0.52918
    with open(path_, "r") as f:
        lines = f.readlines()
    tmp1 = get_line_index("-Fx-", lines)
    tmp2 = get_line_index("Cell Lattice", lines)
    tmp3 = get_line_index("Total stress", lines)

    data = lines[tmp1[0] + 1: tmp2[0] - 1]
    data = np.array([_.split() for _ in data]).astype(float)

    lat = lines[tmp2[0] + 2: tmp3[0] - 1]
    lat = np.array([_.split() for _ in lat]).astype(float)
    cart = data[:, 1:4]

    pos = direct2pos(lat, cart)
    foc = data[:, 4:7] * har2ev / bohr2a

    sts = lines[tmp3[0] + 2: tmp3[0] + 5]
    sts = np.array([_.split() for _ in sts]).astype(float)  # GPa
    vol = np.cross(lat[:, 0], lat[:, 1]) @ lat[:, 2]
    sts = sts * vol * 6.2415093433 * 1e-3  # eV

    ene = get_line_index("Total energy", lines)
    ene = float(lines[ene[0]].split()[3])
    return ene, lat, pos, foc, -sts, ['C'] * len(pos)


if __name__ == '__main__':

    root_dir = argv[1]
    write_dir = argv[2]

    if len(argv) == 3:
        fp_calculator = "vasp"
    else:
        fp_calculator = argv[3]

    if not os.path.exists(write_dir):
        os.mkdir(write_dir)

    if fp_calculator == "vasp":
        car_name = argv[4]
        outcar_dir = find_files(root_dir, car_name)
        for i in outcar_dir:
            xsf_head = get_header(i)
            print(xsf_head)
            try:
                datas = read(i, index=":")
            except Exception as e:
                print(f"{e}")
                continue
            if car_name == "goodStructures_POSCARS":
                datas = read(i, index=":", format="vasp-xdatcar")
            if car_name == "goodStructs_POSCAR":
                datas = read(i, index=":", format="vasp-xdatcar")
            for n, j in enumerate(datas):
                ene, lat, pos, foc, vir, ele = read_from_ase_atoms(j)
                write2my(
                    os.path.join(write_dir, f"{xsf_head}_struct{str(n).zfill(6)}.xsf"),
                    ene_i=ene,
                    lat_i=lat,
                    ele_i=ele,
                    coo_i=pos,
                    foc_i=foc,
                    vir_i=vir)

    if fp_calculator == "qe":
        outfile_dir = find_files(root_dir, "scf.out")
        for i in outfile_dir:
            xsf_head = '_'.join(i.split('.')[-2].split('/'))
            print(xsf_head)
            datas = read(i, index=":", format='espresso-out')
            for n, j in enumerate(datas):
                # print(j)
                ene, lat, pos, foc, vir, ele = read_from_ase_atoms(j)
                write2my(
                    os.path.join(write_dir, f"{xsf_head}_struct{str(n).zfill(6)}.xsf"),
                    ene_i=ene,
                    lat_i=lat,
                    ele_i=ele,
                    coo_i=pos,
                    foc_i=foc,
                    vir_i=vir)

    if fp_calculator == "ares":  # current for carbon only
        car_name = argv[4]
        outcar_dir = find_files(root_dir, car_name)
        for i in outcar_dir:
            xsf_head = get_header(i)
            print(xsf_head)
            datas = [read_ares_scf(i)]

            for n, j in enumerate(datas):
                ene, lat, pos, foc, vir, ele = j
                write2my(
                    os.path.join(write_dir, f"{xsf_head}_struct{str(n).zfill(6)}.xsf"),
                    ene_i=ene,
                    lat_i=lat,
                    ele_i=ele,
                    coo_i=pos,
                    foc_i=foc,
                    vir_i=vir)

