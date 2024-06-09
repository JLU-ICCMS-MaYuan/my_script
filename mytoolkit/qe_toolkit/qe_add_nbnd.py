#!/usr/bin/env python3

import sys
import os
nbnd = sys.argv[1]

def add_or_not(lines):
    '''
    True  表示 需要添加
    False 表示 不需要添加
    '''
    for line in lines:
        if 'nbnd' in line:
            return False
    else:
        return True


with open("scffit.in", 'r') as f2:
    lines2 = f2.readlines()
with open("scf.in", 'r') as f3:
    lines3 = f3.readlines()

if add_or_not(lines2) and add_or_not(lines3):
    os.system(f"sed -i '/&SYSTEM/a  nbnd={nbnd},' scffit.in ")
    os.system(f"sed -i '/&SYSTEM/a  nbnd={nbnd},' scf.in ")
else:
    print("You don't need to add `nbnd` it has existed")
