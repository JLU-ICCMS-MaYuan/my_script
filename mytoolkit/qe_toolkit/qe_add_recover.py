#!/usr/bin/env python3

def add_or_not(lines):
    '''
    True  表示 需要添加
    False 表示 不需要添加
    '''
    for line in lines1:
        if 'recover' in line:
            return False
    else:
        return True
    
with open("split_ph.in", 'r') as f1:
    lines1 = f1.readlines()


if add_or_not(lines1):
    lines1.insert(-1, "  recover=.true.,\n")
    with open("split_ph.in", 'w') as f2:
        f2.writelines(lines1)
else:
    print("You don't need to add `recover=.true.,` it has existed")
