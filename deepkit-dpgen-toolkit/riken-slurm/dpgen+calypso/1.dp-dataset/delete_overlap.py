import os
import numpy as np
from ase import Atoms

def process_directory(directory):
    # 进入目录
    os.chdir(directory)
    
    # 读取文件
    box = np.load('box.npy')         # 形状: (N, 3, 3)
    coord = np.load('coord.npy')     # 形状: (N, M, 3)
    energy = np.load('energy.npy')   # 形状: (N,)
    force = np.load('force.npy')     # 形状: (N, M, 3)
    virial = np.load('virial.npy')   # 形状: (N, 3, 3)

    # N 是快照的数量，M 是每个快照中的原子数量
    N = coord.shape[0]
    M = coord.shape[1]

    threshold_distance = 0.5  # 自定义距离阈值

    # 保留不重叠的快照
    new_box = []
    new_coord = []
    new_energy = []
    new_force = []
    new_virial = []

    for i in range(N):
        # 创建 ASE Atoms 对象
        atoms = Atoms(positions=coord[i].reshape(-1, 3), cell=box[i].reshape(3, 3), pbc=[1, 1, 1])
        distances = atoms.get_all_distances(mic=True)
        
        # 将对角线元素设为1000
        np.fill_diagonal(distances, 1000)
        
        # 检查是否存在小于阈值的距离
        if np.min(distances) > threshold_distance:
            new_box.append(box[i])
            new_coord.append(coord[i])
            new_energy.append(energy[i])
            new_force.append(force[i])
            new_virial.append(virial[i])

    # 将新数据保存为 .npy 文件
    np.save('box.npy', np.array(new_box))
    np.save('coord.npy', np.array(new_coord))
    np.save('energy.npy', np.array(new_energy))
    np.save('force.npy', np.array(new_force))
    np.save('virial.npy', np.array(new_virial))

    print(f"{directory}: 保留了 {len(new_box)} 个快照。")

# 当前路径
current_path = os.getcwd()

# 遍历当前路径下所有目录
for root, dirs, files in os.walk(current_path):
    for directory in dirs:
        if directory.startswith('set'):
            process_directory(os.path.join(root, directory))
