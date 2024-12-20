#!/usr/bin/env python3

#!/usr/bin/env python3
import os
import sys
from ase.io import read, write

print("Please input TRAJ files in the order you wish to merge")
filenames = sys.argv[1:]

atoms_list = []
for filename in filenames:
    # 读取每个traj文件的所有帧
    traj_frames = read(filename, index=':')
    print('filename:{}\nframes found={}\nframes read={}'.format(filename, len(traj_frames), len(traj_frames)))
    atoms_list.extend(traj_frames)

print('Total frames: {}'.format(len(atoms_list)))

# 将所有结构写入一个新的traj文件
write("merged.traj", atoms_list)
