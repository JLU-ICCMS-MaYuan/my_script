#!/usr/bin/env python3

import os
import sys

begin_num = int(sys.argv[1])
end_num   = int(sys.argv[2])

current_path = os.getcwd()
for i in range(begin_num, end_num+1):
    epc_path = os.path.join(current_path, str(i))
    print(f"---------------------------------{i}---------------------------------")
    os.chdir(epc_path)
    os.system(f'ls elph_dir')
    content = os.popen(f"grep freq *dyn{str(i)} | head -n 10").read()
    print(content)
    os.chdir(current_path)
