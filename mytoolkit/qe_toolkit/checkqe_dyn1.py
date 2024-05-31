#!/usr/bin/env python3
import os
import time
import fnmatch

from pathlib import Path

fail_d = []
success_d = []
print("Note: --------------------")
print("    检查qe计算的gamma点的声子是否报错, 如果没报错, 那么是否虚频")
print("    需要你先运行过: checkqe_relax.py. 因为它会直接读取relax-succeeded.dat")

with open("relax-succeeded.dat", "r") as f:
    paths = [ path.strip("\n") for path in f.readlines()]

phonon_failed = []
phonon_succeed = []
phonon_without = []
for path in paths:
    path_1 = Path(path).joinpath(str(1))
    for file in path_1.iterdir():
        if fnmatch.fnmatch(file, '*.dyn1'):
            res = os.system(f"grep -s -m 3 freq {file} ") # 只显示匹配到的文件的指定内容的前三行
            if res != 0: 
                phonon_failed.append(path_1.resolve().__str__())
            else:
                print(file)
                phonon_succeed.append(path_1.resolve().__str__())
            print("检测到dyn1文件，退出"); input()
            break
    else:
        print(path_1)
        print("没有dyn1文件，退出"); input()
        phonon_without.append(path_1.resolve().__str__())

            
with open("phonon-succeeded.dat", "w") as f:
    for line in phonon_succeed:
        f.write(line+"\n")

with open("phonon-failed.dat", "w") as f:
    for line in phonon_failed:
        f.write(line+"\n")

with open("phonon-without.dat", "w") as f:
    for line in phonon_without:
        f.write(line+"\n")