#!/usr/bin/env python3

def delete_or_not(lines):
    '''
    True  表示 需要删除
    False 表示 不需要删除
    '''
    for idx, line in enumerate(lines1):
        if 'recover' in line:
            return idx
    else:
        return None
    
with open("split_ph.in", 'r') as f1:
    lines1 = f1.readlines()


res = delete_or_not(lines1)
if res:
    lines1.pop(res)
    with open("split_ph.in", 'w') as f2:
        f2.writelines(lines1)
    print("Deleted!")
else:
    print("You don't need to delete `recover=.true.,` it doesn't inexist.")
