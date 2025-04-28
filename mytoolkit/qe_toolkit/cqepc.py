#!/usr/bin/env python3

import sys
import os

begin_num = int(sys.argv[1])
end_num   = int(sys.argv[2])

print(begin_num, end_num)

def dyn_done(q):
    for filename in os.listdir(q):
    # 如果文件名包含 "dyn" 关键字
        if   'dyn'+q in filename and os.path.getsize(q+'/'+filename) != 0:
            return True
        elif 'dyn'   in filename and os.path.getsize(q+'/'+filename) != 0:
            return True
    else:
        return False

def imag_freq(freq):
    for f in freq:
        if float(f.split()[4]) < 0:
            return True
    else:
        return False

prefix = os.popen("ls *.dyn0| cut -d '.' -f 1").read().strip()
print(prefix)

print('{:<20} {:<20} {:<20} {:<20}'.format('q_number', 'dyn', 'elph', 'Representation'))
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

    irres_info=''
    if not dyn_done(str(q)):
        dyn_situ = 'inexistence/empty'
        irres_info= os.popen('grep "Representation #" ' +str(q)+ '/split_ph.out | tail -n 1' ).read().strip()
    else:
        if os.path.exists(str(q)+'/'+f'{prefix}.dyn{str(q)}'):
            freq = os.popen('grep freq  '+ str(q)+'/'+f'{prefix}.dyn{str(q)}').readlines()
        elif os.path.exists(str(q)+'/'+f'{prefix}.dyn'):
            freq = os.popen('grep freq  '+ str(q)+'/'+f'{prefix}.dyn').readlines()
        else:
            print('grep freq  '+ str(q)+'/'+f'{prefix}.dyn')
            print("there is no " + str(q)+'/'+f'{prefix}.dyn{str(q)}' + '  or  ' + str(q)+'/'+f'{prefix}.dyn')
        if imag_freq(freq):
            dyn_situ = 'imaginary'
        else:
            dyn_situ = 'done'
    print('{:<20} {:<20} {:<20} {:<}'.format(q, dyn_situ, elph_situ, irres_info))
    situations.append([q, dyn_situ, elph_situ, irres_info])

with open('cqepc.log', 'w') as f:
    f.write('{:<20} {:<20} {:<20} {:<20}\n'.format('q_number', 'dyn', 'elph', 'Representation'))
    for q, dyn_situ, elph_situ, irres_info in situations:
        f.write('{:<20} {:<20} {:<20} {:<}\n'.format(q, dyn_situ, elph_situ, irres_info))
