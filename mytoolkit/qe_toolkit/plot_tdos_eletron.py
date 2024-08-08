#!/usr/bin/env python3
import re
import sys

import matplotlib.pyplot as plt

print("Use it: plot_tdos_eletron.py Nb4H14")
filename = sys.argv[1]

# 常量定义：1 Ry = 13.605693 eV
eV2Ry = 1/13.605693

def get_fermi_energy(file_name):
    with open(file_name, 'r') as file:
        for line in file:
            match = re.search(r'the Fermi energy is\s+(\d+\.\d+)\s*ev', line, re.IGNORECASE)
            if match:
                return float(match.group(1))
    return None

# 从 scffit.out 和 scf.out 文件中提取费米能级
ef_fermi_scffit = get_fermi_energy('scffit.out')
ef_fermi_scf = get_fermi_energy('scf.out')

# 选择非空的费米能级
ef_fermi = ef_fermi_scffit if ef_fermi_scffit is not None else ef_fermi_scf

if ef_fermi is None:
    raise ValueError("无法从 scffit.out 和 scf.out 中提取费米能级")

# 读取 filename.tdos 文件中的数据
energies = []
tdos = []
int_tdos = []

with open(filename+'.tdos', 'r') as file:
    for line in file:
        if not line.startswith('#'):
            data = line.split()
            energies.append(float(data[0]))
            tdos.append(float(data[1]))
            int_tdos.append(float(data[2]))

# 平移能量值，使费米能级位于0
energies_shifted = [e - ef_fermi for e in energies]

# 找到费米能级处的TDOS值
fermi_tdos = None
for i, energy in enumerate(energies_shifted):
    if energy >= 0:
        fermi_tdos = tdos[i]
        print(f"TDOS at fermi = {fermi_tdos} states/eV/f.u.")
        print(f"TDOS at fermi = {fermi_tdos/eV2Ry/2} states/spin/Ry/f.u.")
        break

# 保存平移后的TDOS数据 以 states/eV/Unit Cell (即：states/eV/f.u.) 为单位
with open('Nb4H14_shifted_for_eV.tdos', 'w') as file:
    file.write(f'#  E (eV)   dos(E) (states/eV/f.u.)    Int dos(E)  EFermi = 0.0 eV\n')
    for e, d, i_d in zip(energies_shifted, tdos, int_tdos):
        file.write(f'{e: .3f}  {d: .4E}  {i_d: .4E}\n')

# 保存平移后的TDOS数据 以states/spin/Ry/Unit Cell (即：states/spin/Ry/f.u.) 为单位
with open('Nb4H14_shifted_for_Ry.tdos', 'w') as file:
    file.write(f'#  E (eV)   dos(E) (states/spin/Ry/f.u.)     Int dos(E)  EFermi = 0.0 eV\n')
    for e, d, i_d in zip(energies_shifted, tdos, int_tdos):
        file.write(f'{e: .3f}  {d/eV2Ry/2: .4E}  {i_d/eV2Ry/2: .4E}\n')


# 确定Y轴的上限，增加10%余量以防止标注溢出
y_max = max(tdos) * 1.1

# 绘制数据
plt.figure(figsize=(10, 6))
plt.plot(energies_shifted, tdos, label='TDOS')
plt.axvline(x=0, color='r', linestyle='--', label=f'Fermi level (0 eV)')

# 动态调整标注的位置
if fermi_tdos is not None:
    y_annotation = fermi_tdos + (y_max * 0.05)  # 标注在TDOS值上方
    plt.annotate(f'TDOS at Fermi level: {fermi_tdos}', xy=(0, fermi_tdos), xytext=(0, y_annotation),
                 arrowprops=dict(facecolor='black', shrink=0.05))

plt.xlabel('Energy (eV)')
plt.ylabel('TDOS')
plt.title('Total Density of States')
plt.ylim([0, y_max])  # 设置Y轴的范围
plt.legend()
plt.grid(True)
plt.savefig(filename+'_tdos.png')
plt.show()