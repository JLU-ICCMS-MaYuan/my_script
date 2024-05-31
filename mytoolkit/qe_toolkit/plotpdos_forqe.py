#!/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yyyu200@163.com
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

lw=1.2

F=plt.gcf()
F.clf()
p=plt.subplot(1, 1, 1)

# change by-hand here, read scf input in the future
elem=['La','Y', 'Be', 'H']
ielem=np.array([1,1,2,16],dtype=np.int32) # number of atoms for each element
orb=[['s','p','d','f'],['s','p','d'],['s','p'],['s']]  # projectors for each element
iorb=np.array([6,5,6,16],dtype=np.int32) # number of projectors for each element

num_file=np.dot(ielem,iorb)
nat=np.sum(ielem)
D=[]
N=len(elem)

#scf ATOMIC_POSITIONS should be sorted in the same order as above
count=0
count_at=0
for n in range(N):
    for i in range(ielem[n]):
        for j in range(iorb[n]):
            print(n,i,j,count_at+1,elem[n],j+1,orb[n][j])
            fname='sno.pdos_atm#{}({})_wfc#{}({})'.format(count_at+1,elem[n],j+1,orb[n][j])
            D.append(np.loadtxt(fname,dtype=np.float32))
            count+=1
        count_at+=1

e_fermi=8.7233
line1=plt.plot(D[0][:,0]-e_fermi,D[0][:,1]+D[3][:,1],color='g',linewidth=lw,label='Sn s' )
line2=plt.plot(D[0][:,0]-e_fermi,D[1][:,1]+D[4][:,1],color='r',linewidth=lw,label='Sn p' )
line3=plt.plot(D[0][:,0]-e_fermi,D[2][:,1]+D[5][:,1],color='b',linewidth=lw,label='Sn d' )

line4=plt.plot(D[0][:,0]-e_fermi,D[6][:,1]+D[8][:,1]+D[10][:,1]+D[12][:,1],color='cyan',linewidth=lw,label='O s' )
line5=plt.plot(D[0][:,0]-e_fermi,D[7][:,1]+D[9][:,1]+D[11][:,1]+D[13][:,1],color='grey',linewidth=lw,label='O p' )
plt.xlim([-10,7])
plt.ylim([0,10.5])

plt.ylabel(r'DOS (a.u.)',fontsize=16)
plt.xlabel(r'E (eV) ',fontsize=16)
plt.legend()
plt.savefig('pdos.png',dpi=200)
