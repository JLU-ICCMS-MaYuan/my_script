#!/usr/bin/env python

import numpy as np

with open("hspp.dat", "r") as file:
    lines = file.readlines()

recip_lat = np.array([list(map(float, line.strip().split())) for line in lines[0:3]])
path_name_coords =[ [line.split()[3].strip()] + [list(map(float, line.split()[0:3]))]  for line in lines[3:]]
print("recip_lat: {}".format(recip_lat))
print("path_name_coords: {}".format(path_name_coords))

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

for name, coord in projected_path_name_coords:
    print(f"{name}  {coord}")

    
