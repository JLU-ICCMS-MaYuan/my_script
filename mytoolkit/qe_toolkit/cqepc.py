#!/usr/bin/env python3

import sys
import os

begin_num = int(sys.argv[1])
end_num   = int(sys.argv[2])

print(begin_num, end_num)

def dyn_done(q):
    for filename in os.listdir(q):
    # 如果文件名包含 "dyn" 关键字
        if 'dyn'+q in filename and os.path.getsize(q+'/'+filename) != 0:
            return True
    else:
        return False

def imag_freq(freq):
    for f in freq:
        if float(f.split()[4]) < 0:
            return True
    else:
        return False

print('{:<20} {:<20} {:<20}'.format('q_number', 'dyn', 'elph'))
current_path = os.getcwd()
situations = []
for q in range(begin_num, end_num+1):
    epc_path = os.path.join(current_path, str(q))
    if not os.path.exists(str(q)+'/'+'elph_dir'):
        elph_situ = 'inexistence'
    else:
        if len(os.listdir(str(q)+'/'+'elph_dir')) == 0:
            elph_situ = 'empty'
        else:
            elph_situ = 'done'
    if not dyn_done(str(q)):
        dyn_situ = 'inexistence/empty'
    else:
        freq = os.popen('grep freq  '+ str(q)+'/'+f'*dyn{str(q)}').readlines()
        if imag_freq(freq):
            dyn_situ = 'imaginary'
        else:
            dyn_situ = 'done'
    print('{:<20} {:<20} {:<20}'.format(q, dyn_situ, elph_situ))
    situations.append([q, dyn_situ, elph_situ])

with open('cqepc.log', 'w') as f:
    f.write('{:<20} {:<20} {:<20}\n'.format('q_number', 'dyn', 'elph'))
    for q, dyn_situ, elph_situ in situations:
        f.write('{:<20} {:<20} {:<20}\n'.format(q, dyn_situ, elph_situ))
