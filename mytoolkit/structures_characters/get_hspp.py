#!/usr/bin/env python3
import sys
import argparse

from itertools import chain
from ase.io import read
from ase import Atom
import numpy as np



def interpolate_points(start, end, num_points):
    # return np.linspace(start, end, num_points + 2)[1:-1]  # 不包括端点
    return np.linspace(start, end, num_points + 1)[:-1]  # 包含前一个端点，不包括后一个端点

def write_interpolated_points(path_name_coords, num_points):
    total_points = sum(num_points for _ in range(len(path_name_coords) - 1)) + 1
    with open("path.dat", "w") as f:
        f.write(f"{total_points} crystal\n")  # 添加总点数和“crystal”
        for i in range(len(path_name_coords) - 1):
            start = np.array(path_name_coords[i][1])
            end = np.array(path_name_coords[i + 1][1])
            interpolated = interpolate_points(start, end, num_points)

            # 计算插值点与前一个点之间的距离
            array = (end-start)/num_points
            index = np.nonzero(array)[0]
            # print(array[index])
            distances = np.repeat(np.abs(array[index]), len(interpolated))
            # print(distances)
            
            # 写入插值点及其距离
            for j, point in enumerate(interpolated):
                f.write("{:<10.6f} {:<10.6f} {:<10.6f} {:<10.6f}\n".format(point[0], point[1], point[2], distances[j]))
        f.write("{:<10.6f} {:<10.6f} {:<10.6f} {:<10.6f}\n".format(path_name_coords[-1][1][0], path_name_coords[-1][1][1], path_name_coords[-1][1][2], distances[j]))



def get_hspp(ase_atom:Atom):
    """
    This method is to get high symmetry paths and points
    """ 
    lat       = ase_atom.cell.get_bravais_lattice()
    recip_lat = ase_atom.cell.reciprocal()
    pstring = lat.special_path

    # 获得高对称点路径
    _plist  = [[ p for p in pp if not p.isdigit()] for pp in pstring.split(",")]
    print(f"the high symmetry points path: \n{_plist}")

    print(
        "please input the mode you want, just even input Number like 1 or 2\n",
        "0:  all_points:\n",
        "1:  first_group_points\n",
        "2:  second_group_points\n",
        "n:  n^Th_group_points\n",
        "1 2: first_group_points and second_group_points"
        "...."
        "Nothing to input, directly press ENTER, the default is all_points\n"
        )
    try:
        high_symmetry_type = list(map(int, input().split())) #将输入的整数字符串按照空格划分成列表并分别转化为整数类型并再转化为列表
    except:
        print("what you input is not an integer number, So use the `0:  all_points`")
        high_symmetry_type = [0]

    path_name_list = []
    if "," in pstring:
        if 0 in high_symmetry_type:
            path_name_list = list(chain.from_iterable(_plist))
            print(f"the choosed high symmetry points path is \n {path_name_list}")
        elif 0 not in high_symmetry_type:
            # path_name_list = max(_plist, key=len)
            for hst in high_symmetry_type:
                path_name_list.extend(_plist[hst-1])
            print(f"the choosed high symmetry points path is \n {path_name_list}")
    else:
        path_name_list = [ pp for pp in pstring]

    special_points   = lat.get_special_points()
    path_coords      = [list(special_points[point_name]) for point_name in path_name_list]
    path_name_coords = list(zip(path_name_list, path_coords))


    # 处理高对称点路径
    print("Print Fractional Coordinates of Reciprocal Lattice ! ")
    for name, dirt in path_name_coords:
        print("{:<10.6f} {:<10.6f} {:<10.6f} {:<4}".format(dirt[0], dirt[1], dirt[2], name))
    
    print("The reciprocal lattice (without multiplating `unit_reciprocal_axis`)")
    for vector in recip_lat:
        print("{:<6.3f} {:<6.3f} {:<6.3f} ".format(vector[0], vector[1], vector[2]))

    print("Print projected high symmetry path")
    print("倒格子的单位是 2pi/alat")
    #projected_path_name_coords = [[path_name_coords[0][0], path_name_coords[0][1][0]]]
    projected_path_name_coords = [[path_name_coords[0][0], 0]]
    total_dist = 0
    for idx in range(1, len(path_name_coords)):
        current_name   = path_name_coords[idx][0]
        # current_coords = np.dot(recip_lat, path_name_coords[idx][1])
        # last_coords    = np.dot(recip_lat, path_name_coords[idx-1][1])
        current_coords = np.dot(path_name_coords[idx][1],   recip_lat)
        last_coords    = np.dot(path_name_coords[idx-1][1], recip_lat)
        dist = np.linalg.norm(current_coords-last_coords, 2)
        total_dist += dist
        projected_path_name_coords.append([current_name, total_dist])
    string_names = ' '.join(coord[0] for coord in projected_path_name_coords)
    string_coord = ' '.join(str(np.round(coord[1], 6)) for coord in projected_path_name_coords)
    print(string_names)
    print(string_coord)
    return path_name_coords 

def write4tdep_type(path_name_coords):
    front_path_name_coords = path_name_coords[:-1]
    behind_path_name_coords = path_name_coords[1:]
    path_numbers = len(path_name_coords)-1
    with open("infile.qpoints_dispersion", "w") as f:
        f.write("CUSTOM\n")
        f.write("50\n")
        f.write("{}\n".format(path_numbers))
        for front, behind in zip(front_path_name_coords, behind_path_name_coords):
            f.write('{:<10.8f} {:<10.8f} {:<10.8f}    {:<10.8f} {:<10.8f} {:<10.8f}    {:<5} {:<5}\n'.format(front[1][0], front[1][1], front[1][2], behind[1][0], behind[1][1], behind[1][2], front[0], behind[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process high symmetry paths.')
    parser.add_argument('filename', help='Input filename for ASE atom')
    parser.add_argument('-n', '--num_points', type=int, default=0, help='Number of interpolation points between high symmetry points')
    
    args = parser.parse_args()
    ase_atom = read(args.filename)
    path_name_coords = get_hspp(ase_atom)

    write4tdep_type(path_name_coords)
    if args.num_points > 0:
        write_interpolated_points(path_name_coords, args.num_points)
