import os
import numpy as np
from sys import argv


def read_my_xsf(a_xsf_dir):
    with open(a_xsf_dir, 'r') as file:
        energy = None
        lattice = []

        lines = file.readlines()

        # energy
        for i in lines:
            if 'energy' in i:
                energy = np.array(i.split(' ')[-2], float)
        """print(energy)"""

        # lattice
        for i in lines:
            temp = i.split()
            if len(temp) == 3:
                lattice.append(temp)
        lattice = np.array(lattice, float).reshape(3, 3)
        """print(lattice)"""

        # coordinate & force
        ele_c_f = []
        for i in lines:
            temp = i.split()
            if len(temp) == 7:
                ele_c_f.append(temp)
        ele_c_f = np.array(ele_c_f)
        elements = ele_c_f[:, 0]
        position = ele_c_f[:, 1:4].astype(float)
        force = ele_c_f[:, 4:].astype(float)
        """
        print(elements)
        print(position)
        print(force)
        """
        return energy, lattice, elements, position, force


def get_data(a_xsf_dir):
    energy, lattice, elements, position, force = read_my_xsf(a_xsf_dir)
    types = list()
    for i in elements:
        if i not in types:
            types.append(i)
    num_per_types = []
    coord = []

    for i in types:
        indices_i = [index for index, value in enumerate(elements) if value == i]
        num_per_types.append(len(indices_i))
        coord.append(position[indices_i])

    coord = np.concatenate(coord)
    num_per_types = np.array(num_per_types)

    return lattice, types, num_per_types, coord


def write_POSCAR(file_nam: str,
                 note: str,
                 lattice: np.ndarray,
                 types: list,
                 num_per_types: np.ndarray,
                 coord: np.ndarray,
                 cart: bool = True):
    # check
    assert np.array(lattice).size == 9
    lattice = np.array(lattice).reshape(3, 3)

    assert len(types) == len(num_per_types)

    num_tot = np.sum(num_per_types)
    assert num_tot * 3 == np.array(coord).size
    coord = np.array(coord).reshape(-1, 3)

    f = open(file_nam, "w")
    f.write(note)
    f.write("\n")

    f.write(str(1.0))
    f.write('\n')

    lattice = np.reshape(lattice, [3, 3])
    for row in lattice:
        for col in row:
            f.write(f"{float(col): 22.10f}")
        f.write("\n")

    for i in types:
        f.write(i + "  ")
    f.write('\n')

    for i in num_per_types:
        f.write(str(i) + "  ")
    f.write('\n')

    if cart:
        f.write("Cartesian\n")
    else:
        f.write("Direct\n")

    for row in coord:
        for col in row:
            f.write(f"{float(col): 22.8f}")
        f.write("\n")

    f.close()
    pass


if __name__ == '__main__':
    src_dir = argv[1]
    out_dir = argv[2]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    src_file = [i for i in os.listdir(src_dir) if i[-4:] == ".xsf"]
    fp_name = "fp"

    for file_i in range(len(src_file)):
        src_file_dir = os.path.join(src_dir, src_file[file_i])
        middle_name = fp_name + str(file_i).zfill(4)
        out_file_dir = os.path.join(out_dir, middle_name, "POSCAR")
        if not os.path.exists(os.path.join(out_dir, middle_name)):
            os.mkdir(os.path.join(out_dir, middle_name))

        lattice, types, num_per_types, coord = get_data(src_file_dir)
        write_POSCAR(out_file_dir,
                     "ACNN OOS conf.  " + " ".join(types),
                     lattice,
                     types,
                     num_per_types,
                     coord,
                     cart=True)

