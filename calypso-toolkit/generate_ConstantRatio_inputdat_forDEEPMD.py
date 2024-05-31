#!/usr/bin/env python3
import sys
import os
import copy
import shutil
from itertools import product
from math import pi, pow

def get_volume(rcore:list[float], num:list[int]) -> float:
    v = [4*pow(r*0.529, 3)*pi*n/3 for r, n in zip(rcore, num)]
    return sum(v)*1.5

def get_distance(Atomicradii:list[float]) -> list[list[float], list[float], list[float]]:
    rcore   = [r * 0.529 * 0.7 for r in Atomicradii]
    dists_m = []
    for r, rs in zip(rcore, [rcore]*len(rcore)): # [rcore]*len(rcore) 初步建立距离矩阵
        dists_m.append([rr+r for rr in rs])
    return dists_m

def analyses_inputparas(inputparas: list[str]):
    inputdatas = []
    NameOfAtoms = []
    AtomicNumber = []
    Atomicradii = []
    ctrlRange = []
    for data in inputparas:
        data = data.split()
        assert len(data) == 5
        inputdatas.append(data)
        NameOfAtoms.append(str(data[0]))
        AtomicNumber.append(int(data[1]))
        Atomicradii.append(float(data[2]))
        ctrlRange.append(range(int(data[3]), int(data[4])+1))
    SystemName = '-'.join(NameOfAtoms)
    NumberOfSpecies = len(NameOfAtoms)
    NumberOfAtoms = list(product(*ctrlRange))
    Volumes = [get_volume(Atomicradii, num) for num in NumberOfAtoms]
    DistanceOfIon = get_distance(Atomicradii)
    # print(SystemName, NumberOfSpecies, NameOfAtoms, AtomicNumber, NumberOfAtoms, Volumes, DistanceOfIon)
    return SystemName, NumberOfSpecies, NameOfAtoms, AtomicNumber, NumberOfAtoms, Volumes, DistanceOfIon


def gen_one_inputdat(
        inputdatfile:list[str], 
        SystemName:str,
        NumberOfSpecies:int,
        NameOfAtoms:list[str], 
        AtomicNumber:list[int], 
        NumberOfAtoms:list[int], 
        Volume:float, 
        ):
    new_inputfile = copy.deepcopy(inputdatfile)

    for idx, data in enumerate(new_inputfile):
        if "SystemName" in data:
            new_inputfile[idx] = f"SystemName = {SystemName}\n"
        elif "NumberOfSpecies" in data: # 元素个数
            new_inputfile[idx] = f"NumberOfSpecies = {NumberOfSpecies}\n"
        elif "NameOfAtoms" in data:   # 元素名
            new_inputfile[idx] = f"NameOfAtoms = {' '.join(NameOfAtoms)}\n"
        elif "AtomicNumber" in data:  # 原子序数
            new_inputfile[idx] = f"AtomicNumber = {' '.join(map(str, AtomicNumber))}\n"
        elif "NumberOfAtoms" in data: # 配比
            new_inputfile[idx] = f"NumberOfAtoms = {' '.join(map(str, NumberOfAtoms))}\n"
        elif "Volume" in data:
            new_inputfile[idx] = f"Volume = {Volume}\n"   
        elif "MaxNumAtom" in data:
            new_inputfile[idx] = f"MaxNumAtom = {sum(NumberOfAtoms)*4}\n"  

    return new_inputfile

if __name__ == "__main__":

    info = "Reference command:\npython generate_inputdat.py 'Ce 58 2.55 1 2' 'Sc 21 2.5 1 2' 'H 1 1.1 20 30'"
    print(info)
    print("Ce: AtomicName, 58: AtomicNumber, 2.55: AtomicRadius, 1:min range, 4: max range")

    # 读取当前目录下的input.dat文件
    with open('input.dat') as f:
        inputdatfile = f.readlines()
    
    # 分析输入的参数, 重点是 NumberOfAtoms, 它是
    SystemName,      \
    NumberOfSpecies, \
    NameOfAtoms,     \
    AtomicNumber,    \
    NumberOfAtoms,   \
    Volumes,         \
    DistanceOfIon  = analyses_inputparas(sys.argv[1:])
    
    # 对分析出的每一种化学式建立相应的目录并生成相应的input.dat文件
    serials = 0
    calypso_paths = []
    for NumberOfAtom, Volume in zip(NumberOfAtoms, Volumes):
        serials += 1

        # 准备目录名并且创建目录，并将目录名写入fixed_comp.name文件中，例如：1.fixed_La1Ce3H20
        dirname = ''
        for name, num in zip(NameOfAtoms, NumberOfAtom):
            dirname += name+str(num)
        dirname =  str(serials)+".fixed_"+dirname
        calypso_path = os.path.join(os.path.curdir, dirname)
        if not os.path.exists(calypso_path):
            os.mkdir(calypso_path)
        calypso_paths.append(os.path.basename(calypso_path)+'\n')
        # 获得修改后的input.dat文件 并将其 写入相应的目录中
        new_inputfile = gen_one_inputdat(
            inputdatfile,
            SystemName, 
            NumberOfSpecies, 
            NameOfAtoms, 
            AtomicNumber, 
            NumberOfAtom, 
            Volume, 
            )
        new_inputdat = os.path.join(calypso_path, 'input.dat')
        with open(new_inputdat, 'w') as f:
            f.writelines(new_inputfile)
        
        # 将当前目录下相应的POTCAR, calypso.x, INCAR_*, 
        os.system(f"cp POTCAR    {calypso_path}")
        os.system(f"cp INCAR_*   {calypso_path}")
        os.system(f"cp calypso.x {calypso_path}")
        os.system(f"cp all.pb    {calypso_path}")
        os.system(f"cp calypso_check_outcar.py  {calypso_path}")
        os.system(f"cp calypso_run_opt.py       {calypso_path}")
        os.system(f"cp caly.slurm               {calypso_path}")
        os.system(f"cp submit.sh  {calypso_path}")

    with open("fixed_comp.name", 'w') as f:
        f.writelines(calypso_paths)
