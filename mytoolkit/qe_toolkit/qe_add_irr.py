#!/usr/bin/env python3

import sys

start_irr, last_irr = sys.argv[1:3]

def add_or_not(lines):
    '''
    True  表示 需要添加
    False 表示 不需要添加
    '''
    for line in lines:
        if 'start_irr' in line or 'last_irr' in line:
            return False
    else:
        return True


with open("split_ph.in", 'r') as f2:
    lines2 = f2.readlines()

if add_or_not(lines2):
    lines2.insert(-1, "  start_irr={},\n".format(start_irr))
    lines2.insert(-1, "  last_irr={},\n".format(last_irr))    
    with open("split_ph.in", 'w') as f22:
        f22.writelines(lines2)
    print("Added!")
else:
    print("You don't need to add `start_irr` and `last_irr` it has existed")
