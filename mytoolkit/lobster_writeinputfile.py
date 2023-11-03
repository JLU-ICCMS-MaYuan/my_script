#!/usr/bin/env python
import os 
import sys

from ase.io.pov import get_bondpairs
from ase.io import read

filename = sys.argv[1]
COHPstartEnergy = sys.argv[2]
COHPendEnergy = sys.argv[3]

# 获取成键原子对
atoms = read(filename)
bondpairs_raw = get_bondpairs(atoms)
bondpairs = []
for bp in bondpairs_raw:
    tmp = list(bp[:2])
    tmp.sort()
    if tmp not in bondpairs and tmp[0] != tmp[1]:
        bondpairs.append(tmp)

with open('labels', 'w') as f:
    for bp in bondpairs:
        f.write(f'{atoms.get_chemical_symbols()[bp[0]]}[{bp[0]}]-{atoms.get_chemical_symbols()[bp[1]]}[{bp[1]}]\n')



# 写 lobsterin 输入文件
with open('lobsterin', 'w') as f:
    # 定义能量范围（以费米能级为零点）
    f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
    f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
    f.write('usebasisset pbeVaspFit2015\n') # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    f.write('useRecommendedBasisFunctions\n')
    f.write('gaussianSmearingWidth 0.05\n') # 使用 Gaussian smearing 的话要给出展宽
    # 最后将前一节中获得的成键原子对逐行写入（此处要给原子序号加1，因为python中数组下标是从0开始的
    for bp in bondpairs:
        f.write(f'cohpbetween atom {bp[0]+1} and atom {bp[1]+1}\n')

os.system('$LOBSTER_COMMAND')