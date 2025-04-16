# convert VASP OUTCAT to lammps .data (atomic)
#
# usage:
# python "this file" OUTCAR.input data.output
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ###
import numpy as np
from sys import argv


def anls_poscar(file_name: str):
    with open(file_name, 'r') as f:
        lines = f.readlines()

    bloom_factor = float(lines[1])
    lat = np.array([i.split() for i in lines[2: 5]]).astype(float)
    lat = bloom_factor * lat

    ele_and_num = lines[5][:-1] + ":" + lines[6][:-1]

    pos_d = None
    if (lines[7][0] == 'd') or (lines[7][0] == 'D'):
        # direct coord.
        pos_d = np.array([i.split() for i in lines[8:]]).astype(float)
    elif (lines[7][0] == 'c') or (lines[7][0] == 'C'):
        # cartesian
        pos_c = np.array([i.split() for i in lines[8:]]).astype(float)
        pos_d = direct2pos(lat, pos_c, False)
    return lat, pos_d, ele_and_num


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


def clip_ny_norm(x):
    return x / np.linalg.norm(x)


def cart2lams(lat):
    """
    vasp (POSCAR) Coordinate System:
    O = [[ax, ay, az],
         [bx, by, bz],
         [cx, cy, cz]]


    :return:
    lammps Coordinate System:
    O' = [[xhi - xlo,  0,  0],
          [xy, yhi - ylo,  0],
          [xz, yz, zhi - zlo]]
    """
    ax = (lat[0] ** 2).sum() ** 0.5
    a_hat = clip_ny_norm(lat[0])
    bx = (lat[1] * a_hat).sum()
    by = np.linalg.norm(np.cross(a_hat, lat[1]))
    cx = (lat[2] * a_hat).sum()
    cy = (lat[2] * np.cross(clip_ny_norm(np.cross(lat[0], lat[1])), a_hat)).sum()
    cz = np.abs((lat[2] * clip_ny_norm(np.cross(lat[0], lat[1]))).sum())
    new = np.array([
        [ax, 0, 0],
        [bx, by, 0],
        [cx, cy, cz]
    ])
    return new


def write_lmp_atomic_data(filename, lmp_lat, pos_c, strr):
    new_coo = pos_c
    Z_of_type = strr.split(':')[0].split()
    n_type = strr.split(':')[1].split()
    new_ele = [[i] * int(n_type[i - 1]) for i in range(1, len(n_type) + 1)]
    new_ele = sum(new_ele, [])
    # write
    with open(filename, "w") as f:
        f.writelines("%s (written by ARES-NNP)\n\n" % str(filename + strr))
        f.writelines("%i    atoms\n" % len(new_coo))
        f.writelines("%i    atom types\n" % len(Z_of_type))
        f.writelines("0.0    %.17f  xlo xhi\n" % lmp_lat[0][0])
        f.writelines("0.0    %.17f  ylo yhi\n" % lmp_lat[1][1])
        f.writelines("0.0    %.17f  zlo zhi\n" % lmp_lat[2][2])
        f.writelines(
            "    %18.12f    %18.12f    %18.12f  xy xz yz\n\n\n" % (lmp_lat[1][0], lmp_lat[2][0], lmp_lat[2][1]))
        f.writelines("Atoms\n\n")

        index_list = np.arange(1, len(new_coo) + 1)
        for i, ele, coo in zip(index_list, new_ele, new_coo):
            f.writelines("{0:>6} {1:>3} {2:20.12f} {3:20.12f} {4:20.12f}\n".format(i, ele, coo[0], coo[1], coo[2]))


if __name__ == '__main__':
    filename = argv[2]
    lat, pos, strr = anls_poscar(argv[1])

    crs_para = np.linalg.norm(lat, axis=1)
    exg = np.argsort(crs_para)[::-1]
    lat = lat[exg, :]
    pos = pos[:, exg]

    lmp_lat = cart2lams(lat)
    inv = np.linalg.inv(lat) @ lmp_lat

    pos_c = direct2pos(lat, pos, True) @ inv
    write_lmp_atomic_data(filename, lmp_lat, pos_c, strr)
