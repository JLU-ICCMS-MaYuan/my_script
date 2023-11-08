#!/usr/bin/env python
import os
import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import itertools
knum=50
type=0


def cross(a,b):
    c=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    return c
def dot(a,b):
    c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    return c
def axn(a,b):
    c=[a[0]*b,a[1]*b,a[2]*b]
    return c
def cxa(a,b):
    c=[a[0][0]*b[0]+a[0][1]*b[1]+a[0][2]*b[2],a[1][0]*b[0]+a[1][1]*b[1]+a[1][2]*b[2],a[2][0]*b[0]+a[2][1]*b[1]+a[2][2]*b[2]]
    return c
def mol(a):
    c=(a[0] ** 2+a[1] ** 2+a[2] ** 2) ** 0.5
    return c
def cxc(a,b):
    c=[[],[],[]]
    for i in range(3):
        for j in range(3):
            c[i].append(a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j])
    return c
def printc(a):
    print('%2s %2s %2s\n%2s %2s %2s\n%2s %2s %2s'%(a[0][0],a[0][1],a[0][2],a[1][0],a[1][1],a[1][2],a[2][0],a[2][1],a[2][2]))

os.system('phonopy --symmetry --tolerance=1e-2 POSCAR > sym')
os.system('cp BPOSCAR BPOSCAR0')
nrotate=0
nlattice=1
lattice=[[1,0,0],[0,1,0],[0,0,1]]
rotate=[[1,0,0],[0,1,0],[0,0,1]]

f1 = open('sym','r+')
l1 = f1.read()
allline1 = l1.split("\n")
group = allline1[1].split()[1]
try:
    num_group = int(allline1[1].split()[2][1:-1])
except Exception as e:
    print('phonopy_version: \'1.13.2\'')
    num_group = int(allline1[2].split()[1])


print('groupnum: '+str(num_group))
f1.close()
wn=0
while type==0:



    f = open('BPOSCAR','r+')
    l = f.read()
    allline = l.split("\n")
    name = allline[0]
    ax,ay,az = [ float(x) for x in allline[2].split() ]
    bx,by,bz = [ float(x) for x in allline[3].split() ]
    cx,cy,cz = [ float(x) for x in allline[4].split() ]

    try:
        l_n_atom = [ int(x) for x in allline[5].split() ]
    except Exception as e:
        wn=1    #withname
        l_n_atom = [ int(x) for x in allline[5+wn].split() ]
    # l_atom = allline[0].split()
    n_atom = sum(l_n_atom)
    pos = []
    is_Sd = 0 if len(allline[6+wn].split())==1 else 1
    for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
        pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
        pos.append([pos_a,pos_b,pos_c])
    f.close()
    # n_element = len(l_atom)
    # formula = ''.join([l_atom[i_element]+str(l_n_atom[i_element]) for i_element in range(n_element)])

    a = (ax ** 2 + ay ** 2 + az ** 2) ** 0.5
    b = (bx ** 2 + by ** 2 + bz ** 2) ** 0.5
    c = (cx ** 2 + cy ** 2 + cz ** 2) ** 0.5
    ra = [ax,ay,az]
    rb = [bx,by,bz]
    rc = [cx,cy,cz]
    a=mol(ra)
    b=mol(rb)
    c=mol(rc)
    ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
    kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
    kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
    b1=mol(ka)
    b2=mol(kb)
    b3=mol(kc)

    alpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
    beta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
    gamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
    alpha=math.acos(dot(rb,rc)/mol(rb)/mol(rc))
    beta=math.acos(dot(rc,ra)/mol(rc)/mol(ra))
    gamma=math.acos(dot(ra,rb)/mol(ra)/mol(rb))
    k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
    k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
    k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi

    if num_group<=2:
        system = 'TRI'
        print('system: '+system)
        systemnum=14
        if k_alpha>90 and k_beta>90 and k_gamma>90 and k_gamma==min([k_alpha,k_beta,k_gamma]):
            type = 1
            kpath=[[[0.5,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.5,0.0]],[[0.5,0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.5]],[[0.5,0.0,0.5],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.5,0.5]],[[0.5,0.5,0.5],[0.0,0.0,0.0]]]
            kmark=[['$X$','$\\Gamma$'],['$\\Gamma$','$Y$'],['$L$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$N$','$\\Gamma$'],['$\\Gamma$','$M$'],['$R$','$\\Gamma$']]
        elif k_alpha<90 and k_beta<90 and k_gamma<90 and k_gamma==max([k_alpha,k_beta,k_gamma]):
            type = 2
            kpath=[[[0.0,-0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.5,0.0,0.0]],[[0.5,-0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[-0.5,0.0,0.5]],[[-0.5,-0.5,0.5],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.5]],[[0.0,-0.5,0.5],[0.0,0.0,0.0]]]
            kmark=[['$X$','$\\Gamma$'],['$\\Gamma$','$Y$'],['$L$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$N$','$\\Gamma$'],['$\\Gamma$','$M$'],['$R$','$\\Gamma$']]
        elif k_alpha>90 and k_beta>90 and k_gamma==90:
            type = 3
            kpath=[[[0.5,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.5,0.0]],[[0.5,0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.5]],[[0.5,0.0,0.5],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.5,0.5]],[[0.5,0.5,0.5],[0.0,0.0,0.0]]]
            kmark=[['$X$','$\\Gamma$'],['$\\Gamma$','$Y$'],['$L$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$N$','$\\Gamma$'],['$\\Gamma$','$M$'],['$R$','$\\Gamma$']]
        elif k_alpha<90 and k_beta<90 and k_gamma==90:
            type = 4
            kpath=[[[0.0,-0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.5,0.0,0.0]],[[0.5,-0.5,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[-0.5,0.0,0.5]],[[-0.5,-0.5,0.5],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.5]],[[0.0,-0.5,0.5],[0.0,0.0,0.0]]]
            kmark=[['$X$','$\\Gamma$'],['$\\Gamma$','$Y$'],['$L$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$N$','$\\Gamma$'],['$\\Gamma$','$M$'],['$R$','$\\Gamma$']]

    elif num_group in [5,8,9,12,15]:
        system = 'MCLC'
        print('system: '+system)
        systemnum=13
        tmp=a
        a=b
        b=tmp
        tmp=beta0
        beta0=90
        gamma0=90
        alpha0=180-tmp
        tmp=beta
        beta=math.pi/2
        gamma=math.pi/2
        alpha=math.pi-tmp
        f = open('PPOSCAR','r+')
        l = f.read()
        allline = l.split("\n")
        name = allline[0]
        ax,ay,az = [ float(x) for x in allline[2].split() ]
        bx,by,bz = [ float(x) for x in allline[3].split() ]
        cx,cy,cz = [ float(x) for x in allline[4].split() ]
        l_n_atom = [ int(x) for x in allline[5+wn].split() ]
        # l_atom = allline[0].split()
        n_atom = sum(l_n_atom)
        pos = []
        is_Sd = 0 if len(allline[6+wn].split())==1 else 1
        for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
            pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
            pos.append([pos_a,pos_b,pos_c])
        f.close()
        newpos=[]
        for i_natom in range(n_atom):
            newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])
        w_POSCAR = open('POSCARw','w')
        w_POSCAR.write('%-s\n'%(name))
        w_POSCAR.write('1.0000\n')
        lattice=[[-1,0,0],[0,-1,0],[0,0,1]]
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-ax,-ay,-az))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-bx,-by,-bz))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(cx,cy,cz))
        for l_wPOSCAR in range(5,7+wn+is_Sd):
            w_POSCAR.write(allline[l_wPOSCAR]+'\n')
        for l_wPOSCAR in range(n_atom):
            w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(1-newpos[l_wPOSCAR][0],1-newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
        w_POSCAR.close()
        os.system('mv BPOSCAR oldPOSCAR')
        os.system('mv POSCARw BPOSCAR')
        ra = [-ax,-ay,-az]
        rb = [-bx,-by,-bz]
        rc = [cx,cy,cz]
        pa=mol(ra)
        pb=mol(rb)
        pc=mol(rc)
        palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
        pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
        pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
        ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
        kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
        kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
        b1=mol(ka)
        b2=mol(kb)
        b3=mol(kc)
        k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
        k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
        k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi

        if k_gamma>90:
            type = 1
            xi = (2-b*math.cos(alpha)/c)/(4*math.sin(alpha)**2)
            eta = 1.0/2+2*xi*c*math.cos(alpha)/b
            psi = 3.0/4-a**2/(4*b**2*math.sin(alpha)**2)
            phi = psi+(3.0/4-psi)*b*math.cos(alpha)/c
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$F$'],['$F$','$L$'],['$L$','$I$'],\
                   ['$I_{1}$','$Z$'],['$Z$','$F_{1}$'],['$Y$','$X_{1}$'],\
                   ['$X$','$\\Gamma$'],['$\\Gamma$','$N$'],['$M$','$\\Gamma$']]
            kpath=[[[0,0,0],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,0],[1-xi,1-xi,1-eta]],[[1-xi,1-xi,1-eta],[1.0/2,1.0/2,1.0/2]],[[1.0/2,1.0/2,1.0/2],[phi,1-phi,1.0/2]],\
                   [[1-phi,phi-1,1.0/2],[0,0,1.0/2]],[[0,0,1.0/2],[xi,xi,eta]],[[1.0/2,1.0/2,0],[psi,1-psi,0]],\
                   [[1-psi,psi-1,0],[0,0,0]],[[0,0,0],[1.0/2,0,0]],[[1.0/2,0,1.0/2],[0,0,0]]]
        if k_gamma==90:
            type = 2
            xi = (2-b*math.cos(alpha)/c)/(4*math.sin(alpha)**2)
            eta = 1.0/2+2*xi*c*math.cos(alpha)/b
            psi = 3.0/4-a**2/(4*b**2*math.sin(alpha)**2)
            phi = psi+(3.0/4-psi)*b*math.cos(alpha)/c
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$F$'],['$F$','$L$'],\
                   ['$L$','$I$'],['$I_{1}$','$Z$'],['$Z$','$F_{1}$'],\
                   ['$N$','$\\Gamma$'],['$\\Gamma$','$M$']]
            kpath=[[[0,0,0],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,0],[1-xi,1-xi,1-eta]],[[1-xi,1-xi,1-eta],[1.0/2,1.0/2,1.0/2]],\
                   [[1.0/2,1.0/2,1.0/2],[phi,1-phi,1.0/2]],[[1-phi,phi-1,1.0/2],[0,0,1.0/2]],[[0,0,1.0/2],[xi,xi,eta]],\
                   [[1.0/2,0,0],[0,0,0]],[[0,0,0],[1.0/2,0,1.0/2]]]
        if k_gamma<90 and b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2<=0.999:
            print(b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2)
            type = 3
            mu = (1+b**2/a**2)/4
            delta = b*c*math.cos(alpha)/(2*a**2)
            zeta = mu-1.0/4+(1-b*math.cos(alpha)/c)/(4*math.sin(alpha)**2)
            eta = 1.0/2+2*zeta*c*math.cos(alpha)/b
            phi = 1+zeta-2*mu
            psi = eta-2*delta
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$F$'],['$F$','$H$'],\
                   ['$H$','$Z$'],['$Z$','$I$'],['$I$','$F_{1}$'],\
                   ['$H_{1}$','$Y_{1}$'],['$Y_{1}$','$X$'],['$X$','$\\Gamma$'],\
                   ['$\\Gamma$','$N$'],['$M$','$\\Gamma$']]
            kpath=[[[0,0,0],[mu,mu,delta]],[[mu,mu,delta],[1-phi,1-phi,1-psi]],[[1-phi,1-phi,1-psi],[zeta,zeta,eta]],\
                   [[zeta,zeta,eta],[0,0,1.0/2]],[[0,0,1.0/2],[1.0/2,-1.0/2,1.0/2]],[[1.0/2,-1.0/2,1.0/2],[phi,phi-1,psi]],\
                   [[1-zeta,-zeta,1-eta],[1-mu,-mu,-delta]],[[1-mu,-mu,-delta],[1.0/2,-1.0/2,0]],[[1.0/2,-1.0/2,0],[0,0,0]],\
                   [[0,0,0],[1.0/2,0,0]],[[1.0/2,0,1.0/2],[0,0,0]]]
        if k_gamma<90 and 0.999<b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2<1.001:
            print(b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2)
            type = 4
            mu = (1+b**2/a**2)/4
            delta = b*c*math.cos(alpha)/(2*a**2)
            zeta = mu-1.0/4+(1-b*math.cos(alpha)/c)/(4*math.sin(alpha)**2)
            eta = 1.0/2+2*zeta*c*math.cos(alpha)/b
            phi = 1+zeta-2*mu
            psi = eta-2*delta
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$F$'],['$F$','$H$'],\
                   ['$H$','$Z$'],['$Z$','$I$'],['$H_{1}$','$Y_{1}$'],\
                   ['$Y_{1}$','$X$'],['$X$','$\\Gamma$'],\
                   ['$\\Gamma$','$N$'],['$M$','$\\Gamma$']]
            kpath=[[[0,0,0],[mu,mu,delta]],[[mu,mu,delta],[1-phi,1-phi,1-psi]],[[1-phi,1-phi,1-psi],[zeta,zeta,eta]],\
                   [[zeta,zeta,eta],[0,0,1.0/2]],[[0,0,1.0/2],[1.0/2,-1.0/2,1.0/2]],[[1-zeta,-zeta,1-eta],[1-mu,-mu,-delta]],\
                   [[1-mu,-mu,-delta],[1.0/2,-1.0/2,0]],[[1.0/2,-1.0/2,0],[0,0,0]],\
                   [[0,0,0],[1.0/2,0,0]],[[1.0/2,0,1.0/2],[0,0,0]]]
        if k_gamma<90 and b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2>=1.001:
            print(b*math.cos(alpha)/c+b**2*(math.sin(alpha))**2/a**2)
            type = 5
            zeta = (b**2/a**2+(1-b*math.cos(alpha)/c)/math.sin(alpha)**2)/4
            eta = 1.0/2+2*zeta*c*math.cos(alpha)/b
            mu = eta/2+b**2/(4*a**2)-b*c*math.cos(alpha)/(2*a**2)
            nu = 2*mu-zeta
            rho = 1-zeta*a**2/b**2
            omega = (4*nu-1-b**2*math.sin(alpha)**2/a**2)*c/(2*b*math.cos(alpha))
            delta = zeta*c*math.cos(alpha)/b+omega/2-1.0/4
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$F$'],['$F$','$L$'],\
                   ['$L$','$I$'],['$I_{1}$','$Z$'],['$Z$','$H$'],\
                   ['$H$','$F_{1}$'],['$H_{1}$','$Y_{1}$'],['$Y_{1}$','$X$'],\
                   ['$X$','$\\Gamma$'],['$\\Gamma$','$N$'],['$M$','$\\Gamma$']]
            kpath=[[[0,0,0],[mu,mu,delta]],[[mu,mu,delta],[nu,nu,omega]],[[nu,nu,omega],[1.0/2,1.0/2,1.0/2]],\
                   [[1.0/2,1.0/2,1.0/2],[rho,1-rho,1.0/2]],[[1-rho,rho-1,1.0/2],[0,0,1.0/2]],[[0,0,1.0/2],[zeta,zeta,eta]],\
                   [[zeta,zeta,eta],[1-nu,1-nu,1-omega]],[[1-zeta,-zeta,1-eta],[1-mu,-mu,-delta]],[[1-mu,-mu,-delta],[1.0/2,-1.0/2,0]],\
                   [[1.0/2,-1.0/2,0],[0,0,0]],[[0,0,0],[1.0/2,0,0]],[[1.0/2,0,1.0/2],[0,0,0]]]
        print('pa: %12.8f'%pa)
        print('pb: %12.8f'%pb)
        print('pc: %12.8f\n'%pc)
        print('palpha:  %12.8f'%palpha0)
        print('pbeta:   %12.8f'%pbeta0)
        print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [3,4,6,7,10,11,13,14]:
        system = 'MCL'
        print('system: '+system)
        systemnum=12
        if alpha0<89.9 and b<=c:
            type = 1
            eta = (1-b*math.cos(alpha)/c)/(2*math.sin(alpha)**2)
            nu = 1.0/2-eta*c*math.cos(alpha)/b
            kmark=[['$\\Gamma$','$Y$'],['$Y$','$H$'],['$H$','$C$'],\
                   ['$C$','$E$'],['$E$','$M_{1}$'],['$M_{1}$','$A$'],\
                   ['$A$','$X$'],['$X$','$H_{1}$'],['$M$','$D$'],\
                   ['$D$','$Z$'],['$Y$','$D$']]
            kpath=[[[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[0,eta,1-nu]],[[0,eta,1-nu],[0,1.0/2,1.0/2]],\
                   [[0,1.0/2,1.0/2],[1.0/2,1.0/2,1.0/2]],[[1.0/2,1.0/2,1.0/2],[1.0/2,1-eta,nu]],[[1.0/2,1-eta,nu],[1.0/2,1.0/2,0]],\
                   [[1.0/2,1.0/2,0],[0,1.0/2,0]],[[0,1.0/2,0],[0,1-eta,nu]],[[1.0/2,eta,1-nu],[1.0/2,0,1.0/2]],\
                   [[1.0/2,0,1.0/2],[1.0/2,0,0]],[[0,0,1.0/2],[1.0/2,0,1.0/2]]]
    elif num_group in [146,148,155,160,161,166,167]:
        system = 'RHL'
        print('system: '+system)
        systemnum=11
        f = open('PPOSCAR','r+')
        l = f.read()
        allline = l.split("\n")
        name = allline[0]
        ax,ay,az = [ float(x) for x in allline[2].split() ]
        bx,by,bz = [ float(x) for x in allline[3].split() ]
        cx,cy,cz = [ float(x) for x in allline[4].split() ]
        l_n_atom = [ int(x) for x in allline[5+wn].split() ]
        # l_atom = allline[0].split()
        n_atom = sum(l_n_atom)
        pos = []
        is_Sd = 0 if len(allline[6+wn].split())==1 else 1
        for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
            pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
            pos.append([pos_a,pos_b,pos_c])
        f.close()
        # n_element = len(l_atom)
        # formula = ''.join([l_atom[i_element]+str(l_n_atom[i_element]) for i_element in range(n_element)])

        a = (ax ** 2 + ay ** 2 + az ** 2) ** 0.5
        b = (bx ** 2 + by ** 2 + bz ** 2) ** 0.5
        c = (cx ** 2 + cy ** 2 + cz ** 2) ** 0.5
        ra = [ax,ay,az]
        rb = [bx,by,bz]
        rc = [cx,cy,cz]
        a=mol(ra)
        b=mol(rb)
        c=mol(rc)
        ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
        kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
        kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
        b1=mol(ka)
        b2=mol(kb)
        b3=mol(kc)

        # n_k1 = math.ceil(b1/interval_k)
        # n_k2 = math.ceil(b2/interval_k)
        # n_k3 = math.ceil(b3/interval_k)
        # newpos = []
        #
        # for i_natom in range(n_atom):
        #     newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])

        # #writw KPOINTS
        # w_KPOINTS = open('KPOINTS','w')
        # w_KPOINTS.write('Automatic mesh\n0\nGamma\n  %-3s %-3s %-3s\n  0.  0.  0.\n' %(int(n_k1),int(n_k2),int(n_k3)))
        # w_KPOINTS.close()

        alpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
        beta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
        gamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
        alpha=math.acos(dot(rb,rc)/mol(rb)/mol(rc))
        beta=math.acos(dot(rc,ra)/mol(rc)/mol(ra))
        gamma=math.acos(dot(ra,rb)/mol(ra)/mol(rb))
        k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
        k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
        k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
        os.system('mv BPOSCAR oldPOSCAR')
        os.system('mv PPOSCAR BPOSCAR')
        if alpha0<90:
            type=1
            eta = (1+4*math.cos(alpha))/(2+4*math.cos(alpha))
            nu = 3.0/4-eta/2
            kmark=[['$\\Gamma$','$L$'],['$L$','$B_{1}$'],['$B$','$Z$'],\
                   ['$Z$','$\\Gamma$'],['$\\Gamma$','$X$'],['$Q$','$F$'],\
                   ['$F$','$P_{1}$'],['$P_{1}$','$Z$'],['$L$','$P$']]
            kpath=[[[0,0,0],[1.0/2,0,0]],[[1.0/2,0,0],[1.0/2,1-eta,eta-1]],[[eta,1.0/2,1-eta],[1.0/2,1.0/2,1.0/2]],\
                   [[1.0/2,1.0/2,1.0/2],[0,0,0]],[[0,0,0],[nu,0,-nu]],[[1-nu,nu,0],[1.0/2,1.0/2,0]],\
                   [[1.0/2,1.0/2,0],[1-nu,1-nu,1-eta]],[[1-nu,1-nu,1-eta],[1.0/2,1.0/2,1.0/2]],[[1.0/2,0,0],[eta,nu,nu]]]
        if alpha0>90:
            type=2
            eta = 1.0/(2*math.tan(alpha/2)**2)
            nu = 3.0/4-eta/2
            kmark=[['$\\Gamma$','$P$'],['$P$','$Z$'],['$Z$','$Q$'],\
                   ['$Q$','$\\Gamma$'],['$\\Gamma$','$F$'],['$F$','$P_{1}$'],\
                   ['$P_{1}$','$Q_{1}$'],['$Q_{1}$','$L$'],['$L$','$Z$']]
            kpath=[[[0,0,0],[1-nu,-nu,1-nu]],[[1-nu,-nu,1-nu],[1.0/2,-1.0/2,1.0/2]],[[1.0/2,-1.0/2,1.0/2],[eta,eta,eta]],\
                   [[eta,eta,eta],[0,0,0]],[[0,0,0],[1.0/2,-1.0/2,0]],[[1.0/2,-1.0/2,0],[nu,nu-1,nu-1]],\
                   [[nu,nu-1,nu-1],[1-eta,-eta,-eta]],[[1-eta,-eta,-eta],[1.0/2,0,0]],[[1.0/2,0,0],[1.0/2,-1.0/2,1.0/2]]]
    elif num_group in [143,144,145,147,149,150,151,152,153,154,156,157,158,159,162,163,164,165,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194]:
        system = 'HEX'
        print('system: '+system)
        systemnum=10
        type = 1
        kmark=[['$\\Gamma$','$M$'],['$M$','$K$'],['$K$','$\\Gamma$'],\
               ['$\\Gamma$','$A$'],['$A$','$L$'],['$L$','$H$'],\
               ['$H$','$A$'],['$L$','$M$'],['$K$','$H$']]
        kpath=[[[0,0,0],[1.0/2,0,0]],[[1.0/2,0,0],[1.0/3,1.0/3,0]],[[1.0/3,1.0/3,0],[0,0,0]],\
               [[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[1.0/2,0,1.0/2]],[[1.0/2,0,1.0/2],[1.0/3,1.0/3,1.0/2]],\
               [[1.0/3,1.0/3,1.0/2],[0,0,1.0/2]],[[1.0/2,0,1.0/2],[1.0/2,0,0]],[[1.0/3,1.0/3,0],[1.0/3,1.0/3,1.0/2]]]
    elif num_group in [20,21,35,36,37,38,39,40,41,63,64,65,66,67,68]:
        system = 'ORCC'
        print('system: '+system)
        systemnum=9
        if (nlattice==1 or nlattice==5) and nrotate==0 and a<b:
            type = 1
            zeta = (1+a**2/b**2)/4
            kmark=[['$\\Gamma$','$X$'],['$X$','$S$'],['$S$','$R$'],\
                   ['$R$','$A$'],['$A$','$Z$'],['$Z$','$\\Gamma$'],\
                   ['$\\Gamma$','$Y$'],['$Y$','$X_{1}$'],['$X_{1}$','$A_{1}$'],\
                   ['$A_{1}$','$T$'],['$T$','$Y$'],['$Z$','$T$']]
            kpath=[[[0,0,0],[zeta,zeta,0]],[[zeta,zeta,0],[0,1.0/2,0]],[[0,1.0/2,0],[0,1.0/2,1.0/2]],\
                   [[0,1.0/2,1.0/2],[zeta,zeta,1.0/2]],[[zeta,zeta,1.0/2],[0,0,1.0/2]],[[0,0,1.0/2],[0,0,0]],\
                   [[0,0,0],[-1.0/2,1.0/2,0]],[[-1.0/2,1.0/2,0],[-zeta,1-zeta,0]],[[-zeta,1-zeta,0],[-zeta,1-zeta,1.0/2]],\
                   [[-zeta,1-zeta,1.0/2],[-1.0/2,1.0/2,1.0/2]],[[-1.0/2,1.0/2,1.0/2],[-1.0/2,1.0/2,0]],[[0,0,1.0/2],[-1.0/2,1.0/2,1.0/2]]]
            f = open('PPOSCAR','r+')
            l = f.read()
            allline = l.split("\n")
            name = allline[0]
            ax,ay,az = [ float(x) for x in allline[2].split() ]
            bx,by,bz = [ float(x) for x in allline[3].split() ]
            cx,cy,cz = [ float(x) for x in allline[4].split() ]
            l_n_atom = [ int(x) for x in allline[5+wn].split() ]
            # l_atom = allline[0].split()
            n_atom = sum(l_n_atom)
            pos = []
            is_Sd = 0 if len(allline[6+wn].split())==1 else 1
            for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
                pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
                pos.append([pos_a,pos_b,pos_c])
            f.close()
            newpos=[]
            for i_natom in range(n_atom):
                oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
                oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
                oldcpos=oldrpos*oldr
                newctoold=mat(cxc(rotate,lattice))
                newcpos=oldcpos*newctoold.I
                newr=mat([[a/2,-b/2,0],[a/2,b/2,0],[0,0,c]])
                newrpos=newcpos*newr.I
                print(newrpos)
                newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
            for l_wPOSCAR in range(5,7+wn+is_Sd):
                w_POSCAR.write(allline[l_wPOSCAR]+'\n')
            for l_wPOSCAR in range(n_atom):
                w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
            w_POSCAR.close()
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')
            ra = [newr[0,0],newr[0,1],newr[0,2]]
            rb = [newr[1,0],newr[1,1],newr[1,2]]
            rc = [newr[2,0],newr[2,1],newr[2,2]]
            pa=mol(ra)
            pb=mol(rb)
            pc=mol(rc)
            palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
            pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
            pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
            ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
            kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
            kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
            b1=mol(ka)
            b2=mol(kb)
            b3=mol(kc)
            k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
            k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
            k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
            print('pa: %12.8f'%pa)
            print('pb: %12.8f'%pb)
            print('pc: %12.8f\n'%pc)
            print('palpha:  %12.8f'%palpha0)
            print('pbeta:   %12.8f'%pbeta0)
            print('pgamma:  %12.8f\n'%pgamma0)


            # f = open('PPOSCAR','r+')
            # l = f.read()
            # allline = l.split("\n")
            # name = allline[0]
            # ax,ay,az = [ float(x) for x in allline[2].split() ]
            # bx,by,bz = [ float(x) for x in allline[3].split() ]
            # cx,cy,cz = [ float(x) for x in allline[4].split() ]
            # l_n_atom = [ int(x) for x in allline[5].split() ]
            # l_atom = allline[0].split()
            # n_atom = sum(l_n_atom)
            # pos = []
            # is_Sd = 0 if len(allline[6].split())==1 else 1
            # for l_pos in range(7+is_Sd, 7+n_atom+is_Sd):
            #     pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
            #     pos.append([pos_a,pos_b,pos_c])
            # f.close()
            # n_element = len(l_atom)
            # formula = ''.join([l_atom[i_element]+str(l_n_atom[i_element]) for i_element in range(n_element)])
            #
            # ra = [ax,ay,az]
            # rb = [bx,by,bz]
            # rc = [cx,cy,cz]
            # pa=mol(ra)
            # pb=mol(rb)
            # pc=mol(rc)
            # ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
            # kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
            # kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
            # b1=mol(ka)
            # b2=mol(kb)
            # b3=mol(kc)
            #
            # palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
            # pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
            # pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
            # palpha=math.acos(dot(rb,rc)/mol(rb)/mol(rc))
            # pbeta=math.acos(dot(rc,ra)/mol(rc)/mol(ra))
            # pgamma=math.acos(dot(ra,rb)/mol(ra)/mol(rb))
            # k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
            # k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
            # k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
            # os.system('mv BPOSCAR oldPOSCAR')
            # os.system('mv PPOSCAR BPOSCAR')
            # print('pa: %12.8f'%pa)
            # print('pb: %12.8f'%pb)
            # print('pc: %12.8f\n'%pc)
            # print('palpha:  %12.8f'%palpha0)
            # print('pbeta:   %12.8f'%pbeta0)
            # print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [23,24,44,45,46,71,72,73,74]:
        system = 'ORCI'
        print('system: '+system)
        systemnum=8
        if a<b<c:
            type = 1
            zeta = (1+a**2/c**2)/4
            eta = (1+b**2/c**2)/4
            delta = (b**2-a**2)/(4*c**2)
            mu = (a**2+b**2)/(4*c**2)
            kmark=[['$\\Gamma$','$X$'],['$X$','$L$'],['$L$','$T$'],\
                   ['$T$','$W$'],['$W$','$R$'],['$R$','$X_{1}$'],\
                   ['$X_{1}$','$Z$'],['$Z$','$\\Gamma$'],['$\\Gamma$','$Y$'],\
                   ['$Y$','$S$'],['$S$','$W$'],['$L_{1}$','$Y$'],\
                   ['$Y_{1}$','$Z$']]
            kpath=[[[0,0,0],[-zeta,zeta,zeta]],[[-zeta,zeta,zeta],[-mu,mu,1.0/2-delta]],[[-mu,mu,1.0/2-delta],[0,0,1.0/2]],\
                   [[0,0,1.0/2],[1.0/4,1.0/4,1.0/4]],[[1.0/4,1.0/4,1.0/4],[0,1.0/2,0]],[[0,1.0/2,0],[zeta,1-zeta,-zeta]],\
                   [[zeta,1-zeta,-zeta],[1.0/2,1.0/2,-1.0/2]],[[1.0/2,1.0/2,-1.0/2],[0,0,0]],[[0,0,0],[eta,-eta,eta]],\
                   [[eta,-eta,eta],[1.0/2,0,0]],[[1.0/2,0,0],[1.0/4,1.0/4,1.0/4]],[[mu,-mu,1.0/2+delta],[eta,-eta,eta]],\
                   [[1-eta,eta,-eta],[1.0/2,1.0/2,-1.0/2]]]
            f = open('PPOSCAR','r+')
            l = f.read()
            allline = l.split("\n")
            name = allline[0]
            ax,ay,az = [ float(x) for x in allline[2].split() ]
            bx,by,bz = [ float(x) for x in allline[3].split() ]
            cx,cy,cz = [ float(x) for x in allline[4].split() ]
            l_n_atom = [ int(x) for x in allline[5+wn].split() ]
            # l_atom = allline[0].split()
            n_atom = sum(l_n_atom)
            pos = []
            is_Sd = 0 if len(allline[6+wn].split())==1 else 1
            for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
                pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
                pos.append([pos_a,pos_b,pos_c])
            f.close()
            newpos=[]
            for i_natom in range(n_atom):
                oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
                oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
                oldcpos=oldrpos*oldr
                newctoold=mat(cxc(rotate,lattice))
                newcpos=oldcpos*newctoold.I
                newr=mat([[-a/2,b/2,c/2],[a/2,-b/2,c/2],[a/2,b/2,-c/2]])
                newrpos=newcpos*newr.I
                print(newrpos)
                newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
            for l_wPOSCAR in range(5,7+wn+is_Sd):
                w_POSCAR.write(allline[l_wPOSCAR]+'\n')
            for l_wPOSCAR in range(n_atom):
                w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
            w_POSCAR.close()
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')
            ra = [newr[0,0],newr[0,1],newr[0,2]]
            rb = [newr[1,0],newr[1,1],newr[1,2]]
            rc = [newr[2,0],newr[2,1],newr[2,2]]
            pa=mol(ra)
            pb=mol(rb)
            pc=mol(rc)
            palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
            pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
            pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
            ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
            kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
            kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
            b1=mol(ka)
            b2=mol(kb)
            b3=mol(kc)
            k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
            k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
            k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
            print('pa: %12.8f'%pa)
            print('pb: %12.8f'%pb)
            print('pc: %12.8f\n'%pc)
            print('palpha:  %12.8f'%palpha0)
            print('pbeta:   %12.8f'%pbeta0)
            print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [22,42,43,69,70]:
        system = 'ORCF'
        print('system: '+system)
        systemnum=7
        if a<b<c:
            if (1/a**2)*0.999 >= 1/b**2+1/c**2:
                print(1/a**2)
                print(1/b**2+1/c**2)
                type=1
                zeta = (1+a**2/b**2-a**2/c**2)/4
                eta = (1+a**2/b**2+a**2/c**2)/4
                kmark=[['$\\Gamma$','$Y$'],['$Y$','$T$'],['$T$','$Z$'],\
                       ['$Z$','$\\Gamma$'],['$\\Gamma$','$X$'],['$X$','$A_{1}$'],\
                       ['$A_{1}$','$Y$'],['$T$','$X_{1}$'],['$X$','$A$'],\
                       ['$A$','$Z$'],['$L$','$\\Gamma$']]
                kpath=[[[0,0,0],[1.0/2,0,1.0/2]],[[1.0/2,0,1.0/2],[1,1.0/2,1.0/2]],[[1,1.0/2,1.0/2],[1.0/2,1.0/2,0]],\
                       [[1.0/2,1.0/2,0],[0,0,0]],[[0,0,0],[0,eta,eta]],[[0,eta,eta],[1.0/2,1.0/2-zeta,1-zeta]],\
                       [[1.0/2,1.0/2-zeta,1-zeta],[1.0/2,0,1.0/2]],[[1,1.0/2,1.0/2],[1,1-eta,1-eta]],[[0,eta,eta],[1.0/2,1.0/2+zeta,zeta]],\
                       [[1.0/2,1.0/2+zeta,zeta],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,1.0/2],[0,0,0]]]
            if (1/a**2)*1.001 <= 1/b**2+1/c**2:
                type=2
                print(1/a**2)
                print(1/b**2+1/c**2)
                eta = (1+a**2/b**2-a**2/c**2)/4
                delta = (1+b**2/a**2-b**2/c**2)/4
                phi = (1+c**2/b**2-c**2/a**2)/4
                kmark=[['$\\Gamma$','$Y$'],['$Y$','$C$'],['$C$','$D$'],\
                       ['$D$','$X$'],['$X$','$\\Gamma$'],['$\\Gamma$','$Z$'],\
                       ['$Z$','$D_{1}$'],['$D_{1}$','$H$'],['$H$','$C$'],\
                       ['$C_{1}$','$Z$'],['$X$','$H_{1}$'],['$H$','$Y$'],\
                       ['$L$','$\\Gamma$']]
                kpath=[[[0,0,0],[1.0/2,0,1.0/2]],[[1.0/2,0,1.0/2],[1.0/2,1.0/2-eta,1-eta]],[[1.0/2,1.0/2-eta,1-eta],[1.0/2-delta,1.0/2,1-delta]],\
                       [[1.0/2-delta,1.0/2,1-delta],[0,1.0/2,1.0/2]],[[0,1.0/2,1.0/2],[0,0,0]],[[0,0,0],[1.0/2,1.0/2,0]],\
                       [[1.0/2,1.0/2,0],[1.0/2+delta,1.0/2,delta]],[[1.0/2+delta,1.0/2,delta],[1-phi,1.0/2-phi,1.0/2]],[[1-phi,1.0/2-phi,1.0/2],[1.0/2,1.0/2-eta,1-eta]],\
                       [[1.0/2,1.0/2+eta,eta],[1.0/2,1.0/2,0]],[[0,1.0/2,1.0/2],[phi,1.0/2+phi,1.0/2]],[[1-phi,1.0/2-phi,1.0/2],[1.0/2,0,1.0/2]],\
                       [[1.0/2,1.0/2,1.0/2],[0,0,0]]]
            if (1/a**2)*1.001 > 1/b**2+1/c**2 > (1/a**2)*0.999:
                type=3
                print(1/a**2)
                print(1/b**2+1/c**2)
                zeta = (1+a**2/b**2-a**2/c**2)/4
                eta = (1+a**2/b**2+a**2/c**2)/4
                kmark=[['$\\Gamma$','$Y$'],['$Y$','$T$'],['$T$','$Z$'],\
                       ['$Z$','$\\Gamma$'],['$\\Gamma$','$X$'],['$X$','$A_{1}$'],\
                       ['$A_{1}$','$Y$'],['$X$','$A$'],\
                       ['$A$','$Z$'],['$L$','$\\Gamma$']]
                kpath=[[[0,0,0],[1.0/2,0,1.0/2]],[[1.0/2,0,1.0/2],[1,1.0/2,1.0/2]],[[1,1.0/2,1.0/2],[1.0/2,1.0/2,0]],\
                       [[1.0/2,1.0/2,0],[0,0,0]],[[0,0,0],[0,eta,eta]],[[0,eta,eta],[1.0/2,1.0/2-zeta,1-zeta]],\
                       [[1.0/2,1.0/2-zeta,1-zeta],[1.0/2,0,1.0/2]],[[0,eta,eta],[1.0/2,1.0/2+zeta,zeta]],\
                       [[1.0/2,1.0/2+zeta,zeta],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,1.0/2],[0,0,0]]]
            f = open('PPOSCAR','r+')
            l = f.read()
            allline = l.split("\n")
            name = allline[0]
            ax,ay,az = [ float(x) for x in allline[2].split() ]
            bx,by,bz = [ float(x) for x in allline[3].split() ]
            cx,cy,cz = [ float(x) for x in allline[4].split() ]
            l_n_atom = [ int(x) for x in allline[5+wn].split() ]
            # l_atom = allline[0].split()
            n_atom = sum(l_n_atom)
            pos = []
            is_Sd = 0 if len(allline[6+wn].split())==1 else 1
            for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
                pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
                pos.append([pos_a,pos_b,pos_c])
            f.close()
            newpos=[]
            for i_natom in range(n_atom):
                oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
                oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
                oldcpos=oldrpos*oldr
                newctoold=mat(cxc(rotate,lattice))
                newcpos=oldcpos*newctoold.I
                newr=mat([[0,b/2,c/2],[a/2,0,c/2],[a/2,b/2,0]])
                newrpos=newcpos*newr.I
                print(newrpos)
                newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
            for l_wPOSCAR in range(5,7+wn+is_Sd):
                w_POSCAR.write(allline[l_wPOSCAR]+'\n')
            for l_wPOSCAR in range(n_atom):
                w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
            w_POSCAR.close()
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')
            ra = [newr[0,0],newr[0,1],newr[0,2]]
            rb = [newr[1,0],newr[1,1],newr[1,2]]
            rc = [newr[2,0],newr[2,1],newr[2,2]]
            pa=mol(ra)
            pb=mol(rb)
            pc=mol(rc)
            palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
            pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
            pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
            ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
            kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
            kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
            b1=mol(ka)
            b2=mol(kb)
            b3=mol(kc)
            k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
            k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
            k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
            print('pa: %12.8f'%pa)
            print('pb: %12.8f'%pb)
            print('pc: %12.8f\n'%pc)
            print('palpha:  %12.8f'%palpha0)
            print('pbeta:   %12.8f'%pbeta0)
            print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [16,17,18,19,25,26,27,28,29,30,31,32,33,34,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62]:
        system = 'ORC'
        print('system: '+system)
        systemnum=6
        if a<b<c:
            type = 1
            kmark=[['$\\Gamma$','$X$'],['$X$','$S$'],['$S$','$Y$'],\
                   ['$Y$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$Z$','$U$'],\
                   ['$U$','$R$'],['$R$','$T$'],['$T$','$Z$'],\
                   ['$Y$','$T$'],['$U$','$X$'],['$S$','$R$']]
            kpath=[[[0,0,0],[1.0/2,0,0]],[[1.0/2,0,0],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,0],[0,1.0/2,0]],\
                   [[0,1.0/2,0],[0,0,0]],[[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[1.0/2,0,1.0/2]],\
                   [[1.0/2,0,1.0/2],[1.0/2,1.0/2,1.0/2]],[[1.0/2,1.0/2,1.0/2],[0,1.0/2,1.0/2]],[[0,1.0/2,1.0/2],[0,0,1.0/2]],\
                   [[0,1.0/2,0],[0,1.0/2,1.0/2]],[[1.0/2,0,1.0/2],[1.0/2,0,0]],[[1.0/2,1.0/2,0],[1.0/2,1.0/2,1.0/2]]]
    elif num_group in [79,80,82,87,88,97,98,107,108,109,110,119,120,121,122,139,140,141,142]:
        system = 'BCT'
        print('system: '+system)
        systemnum=5
        if a==b:
            if c<a:
                type=1
                eta = (1+c**2/a**2)/4
                kmark=[['$\\Gamma$','$X$'],['$X$','$M$'],['$M$','$\\Gamma$'],\
                       ['$\\Gamma$','$Z$'],['$Z$','$P$'],['$P$','$N$'],\
                       ['$N$','$Z_{1}$'],['$Z_{1}$','$M$'],['$X$','$P$']]
                kpath=[[[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[-1.0/2,1.0/2,1.0/2]],[[-1.0/2,1.0/2,1.0/2],[0,0,0]],\
                       [[0,0,0],[eta,eta,-eta]],[[eta,eta,-eta],[1.0/4,1.0/4,1.0/4]],[[1.0/4,1.0/4,1.0/4],[0,1.0/2,0]],\
                       [[0,1.0/2,0],[-eta,1-eta,eta]],[[-eta,1-eta,eta],[-1.0/2,1.0/2,1.0/2]],[[0,0,1.0/2],[1.0/4,1.0/4,1.0/4]]]
            if c>a:
                type=2
                eta = (1+a**2/c**2)/4
                zeta = a**2/(2*c**2)
                kmark=[['$\\Gamma$','$X$'],['$X$','$Y$'],['$Y$','$\\Sigma$'],\
                       ['$\\Sigma$','$\\Gamma$'],['$\\Gamma$','$Z$'],['$Z$','$\\Sigma_{1}$'],\
                       ['$\\Sigma_{1}$','$N$'],['$N$','$P$'],['$P$','$Y_{1}$'],\
                       ['$Y_{1}$','$Z$'],['$X$','$P$']]
                kpath=[[[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[-zeta,zeta,1.0/2]],[[-zeta,zeta,1.0/2],[-eta,eta,eta]],\
                       [[-eta,eta,eta],[0,0,0]],[[0,0,0],[1.0/2,1.0/2,-1.0/2]],[[1.0/2,1.0/2,-1.0/2],[eta,1-eta,-eta]],\
                       [[eta,1-eta,-eta],[0,1.0/2,0]],[[0,1.0/2,0],[1.0/4,1.0/4,1.0/4]],[[1.0/4,1.0/4,1.0/4],[1.0/2,1.0/2,-zeta]],\
                       [[1.0/2,1.0/2,-zeta],[1.0/2,1.0/2,-1.0/2]],[[0,0,1.0/2],[1.0/4,1.0/4,1.0/4]]]
            f = open('PPOSCAR','r+')
            l = f.read()
            allline = l.split("\n")
            name = allline[0]
            ax,ay,az = [ float(x) for x in allline[2].split() ]
            bx,by,bz = [ float(x) for x in allline[3].split() ]
            cx,cy,cz = [ float(x) for x in allline[4].split() ]
            l_n_atom = [ int(x) for x in allline[5+wn].split() ]
            # l_atom = allline[0].split()
            n_atom = sum(l_n_atom)
            pos = []
            is_Sd = 0 if len(allline[6+wn].split())==1 else 1
            for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
                pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
                pos.append([pos_a,pos_b,pos_c])
            f.close()
            newpos=[]
            for i_natom in range(n_atom):
                oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
                oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
                oldcpos=oldrpos*oldr
                newctoold=mat(cxc(rotate,lattice))
                newcpos=oldcpos*newctoold.I
                newr=mat([[-a/2,b/2,c/2],[a/2,-b/2,c/2],[a/2,b/2,-c/2]])
                newrpos=newcpos*newr.I
                print(newrpos)
                newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
            for l_wPOSCAR in range(5,7+wn+is_Sd):
                w_POSCAR.write(allline[l_wPOSCAR]+'\n')
            for l_wPOSCAR in range(n_atom):
                w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
            w_POSCAR.close()
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')
            ra = [newr[0,0],newr[0,1],newr[0,2]]
            rb = [newr[1,0],newr[1,1],newr[1,2]]
            rc = [newr[2,0],newr[2,1],newr[2,2]]
            pa=mol(ra)
            pb=mol(rb)
            pc=mol(rc)
            palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
            pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
            pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
            ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
            kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
            kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
            b1=mol(ka)
            b2=mol(kb)
            b3=mol(kc)
            k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
            k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
            k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
            print('pa: %12.8f'%pa)
            print('pb: %12.8f'%pb)
            print('pc: %12.8f\n'%pc)
            print('palpha:  %12.8f'%palpha0)
            print('pbeta:   %12.8f'%pbeta0)
            print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [75,76,77,78,81,83,84,85,86,89,90,91,92,93,94,95,96,99,100,101,102,103,104,105,106,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138]:
        system = 'TET'
        print('system: '+system)
        systemnum=4
        if a==b:
            type = 1
            kmark=[['$\\Gamma$','$X$'],['$X$','$M$'],['$M$','$\\Gamma$'],\
                   ['$\\Gamma$','$Z$'],['$Z$','$R$'],['$R$','$A$'],\
                   ['$A$','$Z$'],['$X$','$R$'],['$M$','$A$']]
            kpath=[[[0,0,0],[0,1.0/2,0]],[[0,1.0/2,0],[1.0/2,1.0/2,0]],[[1.0/2,1.0/2,0],[0,0,0]],\
                   [[0,0,0],[0,0,1.0/2]],[[0,0,1.0/2],[0,1.0/2,1.0/2]],[[0,1.0/2,1.0/2],[1.0/2,1.0/2,1.0/2]],\
                   [[1.0/2,1.0/2,1.0/2],[0,0,1.0/2]],[[0,1.0/2,0],[0,1.0/2,1.0/2]],[[1.0/2,1.0/2,0],[1.0/2,1.0/2,1.0/2]]]
    elif num_group in [197,199,204,206,211,214,217,220,229,230]:
        system = 'BCC'
        print('system: '+system)
        systemnum=3
        type = 1
        kmark=[['$\\Gamma$','$H$'],['$H$','$N$'],\
               ['$N$','$\\Gamma$'],['$\\Gamma$','$P$'],\
               ['$P$','$H$'],['$P$','$N$']]
        kpath=[[[0,0,0],[1.0/2,-1.0/2,1.0/2]],[[1.0/2,-1.0/2,1.0/2],[0,0,1.0/2]],\
               [[0,0,1.0/2],[0,0,0]],[[0,0,0],[1.0/4,1.0/4,1.0/4]],\
               [[1.0/4,1.0/4,1.0/4],[1.0/2,-1.0/2,1.0/2]],[[1.0/4,1.0/4,1.0/4],[0,0,1.0/2]]]
        f = open('PPOSCAR','r+')
        l = f.read()
        allline = l.split("\n")
        name = allline[0]
        ax,ay,az = [ float(x) for x in allline[2].split() ]
        bx,by,bz = [ float(x) for x in allline[3].split() ]
        cx,cy,cz = [ float(x) for x in allline[4].split() ]
        l_n_atom = [ int(x) for x in allline[5+wn].split() ]
        # l_atom = allline[0].split()
        n_atom = sum(l_n_atom)
        pos = []
        is_Sd = 0 if len(allline[6+wn].split())==1 else 1
        for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
            pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
            pos.append([pos_a,pos_b,pos_c])
        f.close()
        newpos=[]
        for i_natom in range(n_atom):
            oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
            oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
            oldcpos=oldrpos*oldr
            newctoold=mat(cxc(rotate,lattice))
            newcpos=oldcpos*newctoold.I
            newr=mat([[-a/2,b/2,c/2],[a/2,-b/2,c/2],[a/2,b/2,-c/2]])
            newrpos=newcpos*newr.I
            print(newrpos)
            newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
        w_POSCAR = open('POSCARw','w')
        w_POSCAR.write('%-s\n'%(name))
        w_POSCAR.write('1.0000\n')
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
        for l_wPOSCAR in range(5,7+wn+is_Sd):
            w_POSCAR.write(allline[l_wPOSCAR]+'\n')
        for l_wPOSCAR in range(n_atom):
            w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
        w_POSCAR.close()
        os.system('mv BPOSCAR oldPOSCAR')
        os.system('mv POSCARw BPOSCAR')
        ra = [newr[0,0],newr[0,1],newr[0,2]]
        rb = [newr[1,0],newr[1,1],newr[1,2]]
        rc = [newr[2,0],newr[2,1],newr[2,2]]
        pa=mol(ra)
        pb=mol(rb)
        pc=mol(rc)
        palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
        pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
        pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
        ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
        kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
        kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
        b1=mol(ka)
        b2=mol(kb)
        b3=mol(kc)
        k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
        k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
        k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
        print('pa: %12.8f'%pa)
        print('pb: %12.8f'%pb)
        print('pc: %12.8f\n'%pc)
        print('palpha:  %12.8f'%palpha0)
        print('pbeta:   %12.8f'%pbeta0)
        print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [196,202,203,209,210,216,219,225,226,227,228]:
        system = 'FCC'
        print('system: '+system)
        systemnum=2
        type = 1
        kmark=[['$\\Gamma$','$X$'],['$X$','$W$'],['$W$','$K$'],\
               ['$K$','$\\Gamma$'],['$\\Gamma$','$L$'],['$L$','$U$'],\
               ['$U$','$W$'],['$W$','$L$'],\
               ['$L$','$K$'],['$U$','$X$']]
        kpath=[[[0,0,0],[1.0/2,0,1.0/2]],[[1.0/2,0,1.0/2],[1.0/2,1.0/4,3.0/4]],[[1.0/2,1.0/4,3.0/4],[3.0/8,3.0/8,3.0/4]],\
               [[3.0/8,3.0/8,3.0/4],[0,0,0]],[[0,0,0],[1.0/2,1.0/2,1.0/2]],[[1.0/2,1.0/2,1.0/2],[5.0/8,1.0/4,5.0/8]],\
               [[5.0/8,1.0/4,5.0/8],[1.0/2,1.0/4,3.0/4]],[[1.0/2,1.0/4,3.0/4],[1.0/2,1.0/2,1.0/2]],\
               [[1.0/2,1.0/2,1.0/2],[3.0/8,3.0/8,3.0/4]],[[5.0/8,1.0/4,5.0/8],[1.0/2,0,1.0/2]]]
        f = open('PPOSCAR','r+')
        l = f.read()
        allline = l.split("\n")
        name = allline[0]
        ax,ay,az = [ float(x) for x in allline[2].split() ]
        bx,by,bz = [ float(x) for x in allline[3].split() ]
        cx,cy,cz = [ float(x) for x in allline[4].split() ]
        l_n_atom = [ int(x) for x in allline[5+wn].split() ]
        # l_atom = allline[0].split()
        n_atom = sum(l_n_atom)
        pos = []
        is_Sd = 0 if len(allline[6+wn].split())==1 else 1
        for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
            pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
            pos.append([pos_a,pos_b,pos_c])
        f.close()
        newpos=[]
        for i_natom in range(n_atom):
            oldrpos=mat([[pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]]])
            oldr=mat([[ax,ay,az],[bx,by,bz],[cx,cy,cz]])
            oldcpos=oldrpos*oldr
            newctoold=mat(cxc(rotate,lattice))
            newcpos=oldcpos*newctoold.I
            newr=mat([[0,b/2,c/2],[a/2,0,c/2],[a/2,b/2,0]])
            newrpos=newcpos*newr.I
            print(newrpos)
            newpos.append([newrpos[0,0]-math.floor(newrpos[0,0]),newrpos[0,1]-math.floor(newrpos[0,1]),newrpos[0,2]-math.floor(newrpos[0,2])])
        w_POSCAR = open('POSCARw','w')
        w_POSCAR.write('%-s\n'%(name))
        w_POSCAR.write('1.0000\n')
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[0,0],newr[0,1],newr[0,2]))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[1,0],newr[1,1],newr[1,2]))
        w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(newr[2,0],newr[2,1],newr[2,2]))
        for l_wPOSCAR in range(5,7+wn+is_Sd):
            w_POSCAR.write(allline[l_wPOSCAR]+'\n')
        for l_wPOSCAR in range(n_atom):
            w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
        w_POSCAR.close()
        os.system('mv BPOSCAR oldPOSCAR')
        os.system('mv POSCARw BPOSCAR')
        ra = [newr[0,0],newr[0,1],newr[0,2]]
        rb = [newr[1,0],newr[1,1],newr[1,2]]
        rc = [newr[2,0],newr[2,1],newr[2,2]]
        pa=mol(ra)
        pb=mol(rb)
        pc=mol(rc)
        palpha0=math.acos(dot(rb,rc)/mol(rb)/mol(rc))*180/math.pi
        pbeta0=math.acos(dot(rc,ra)/mol(rc)/mol(ra))*180/math.pi
        pgamma0=math.acos(dot(ra,rb)/mol(ra)/mol(rb))*180/math.pi
        ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
        kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
        kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
        b1=mol(ka)
        b2=mol(kb)
        b3=mol(kc)
        k_alpha=math.acos(dot(kb,kc)/mol(kb)/mol(kc))*180/math.pi
        k_beta=math.acos(dot(kc,ka)/mol(kc)/mol(ka))*180/math.pi
        k_gamma=math.acos(dot(ka,kb)/mol(ka)/mol(kb))*180/math.pi
        print('pa: %12.8f'%pa)
        print('pb: %12.8f'%pb)
        print('pc: %12.8f\n'%pc)
        print('palpha:  %12.8f'%palpha0)
        print('pbeta:   %12.8f'%pbeta0)
        print('pgamma:  %12.8f\n'%pgamma0)
    elif num_group in [195,198,200,201,205,207,208,212,213,215,218,221,222,223,224]:
        system = 'CUB'
        print('system: '+system)
        systemnum=1
        type = 1
        kmark=[['$\\Gamma$','$X$'],['$X$','$M$'],\
               ['$M$','$\\Gamma$'],['$\\Gamma$','$R$'],\
               ['$R$','$X$'],['$M$','$R$']]
        kpath=[[[0,0,0],[0,1.0/2,0]],[[0,1.0/2,0],[1.0/2,1.0/2,0]],\
               [[1.0/2,1.0/2,0],[0,0,0]],[[0,0,0],[1.0/2,1.0/2,1.0/2]],\
               [[1.0/2,1.0/2,1.0/2],[0,1.0/2,0]],[[1.0/2,1.0/2,0],[1.0/2,1.0/2,1.0/2]]]

    print('lattice: '+str(nlattice))
    printc(lattice)
    print('rotate: '+str(nrotate))
    printc(rotate)
    print('matrix to BPOSCAR0:')
    printc(cxc(rotate,lattice))
    print('a: %12.8f'%a)
    print('b: %12.8f'%b)
    print('c: %12.8f\n'%c)
    print('alpha:  %12.8f'%alpha0)
    print('beta:   %12.8f'%beta0)
    print('gamma:  %12.8f\n'%gamma0)
    print('k1: %12.8f'%b1)
    print('k2: %12.8f'%b2)
    print('k3: %12.8f\n'%b3)
    print('kalpha: %12.8f'%k_alpha)
    print('kbeta:  %12.8f'%k_beta)
    print('kgamma: %12.8f'%k_gamma)
    print('type: '+str(type)+'\n-----------------------------------')

    if type==0:
        nrotate+=1
        if nrotate==3:
            nrotate=0
            rotate=[[1,0,0],[0,1,0],[0,0,1]]
            f = open('BPOSCAR0','r+')
            l = f.read()
            allline = l.split("\n")
            name = allline[0]
            ax,ay,az = [ float(x) for x in allline[2].split() ]
            bx,by,bz = [ float(x) for x in allline[3].split() ]
            cx,cy,cz = [ float(x) for x in allline[4].split() ]
            l_n_atom = [ int(x) for x in allline[5+wn].split() ]
            # l_atom = allline[0].split()
            n_atom = sum(l_n_atom)
            pos = []
            is_Sd = 0 if len(allline[6+wn].split())==1 else 1
            for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
                pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
                pos.append([pos_a,pos_b,pos_c])
            f.close()
            newpos=[]
            for i_natom in range(n_atom):
                newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            if nlattice==1:
                lattice=[[1,0,0],[0,-1,0],[0,0,-1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(ax,ay,az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-bx,-by,-bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-cx,-cy,-cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][0],1-newpos[l_wPOSCAR][1],1-newpos[l_wPOSCAR][2]))
            if nlattice==2:
                lattice=[[-1,0,0],[0,1,0],[0,0,-1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-ax,-ay,-az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(bx,by,bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-cx,-cy,-cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(1-newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1],1-newpos[l_wPOSCAR][2]))
            if nlattice==3:
                lattice=[[-1,0,0],[0,-1,0],[0,0,1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-ax,-ay,-az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-bx,-by,-bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(cx,cy,cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(1-newpos[l_wPOSCAR][0],1-newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][2]))
            if nlattice==4:
                lattice=[[0,1,0],[-1,0,0],[0,0,1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(bx,by,bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-ax,-ay,-az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(cx,cy,cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][1],1-newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][2]))
            if nlattice==5:
                lattice=[[0,-1,0],[1,0,0],[0,0,1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-bx,-by,-bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(ax,ay,az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(cx,cy,cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(1-newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][2]))
            if nlattice==6:
                lattice=[[0,1,0],[1,0,0],[0,0,-1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(bx,by,bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(ax,ay,az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-cx,-cy,-cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][1],newpos[l_wPOSCAR][0],1-newpos[l_wPOSCAR][2]))
            if nlattice==7:
                lattice=[[0,-1,0],[-1,0,0],[0,0,-1]]
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-bx,-by,-bz))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-ax,-ay,-az))
                w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(-cx,-cy,-cz))
                for l_wPOSCAR in range(5,7+wn+is_Sd):
                    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
                for l_wPOSCAR in range(n_atom):
                    w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(1-newpos[l_wPOSCAR][1],1-newpos[l_wPOSCAR][0],1-newpos[l_wPOSCAR][2]))
            w_POSCAR.close()
            nlattice+=1
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')


        else:
            rt=[[0,0,1],[1,0,0],[0,1,0]]
            rotate=cxc(rt,rotate)
            newpos=[]
            for i_natom in range(n_atom):
                newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])
            w_POSCAR = open('POSCARw','w')
            w_POSCAR.write('%-s\n'%(name))
            w_POSCAR.write('1.0000\n')
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(cx,cy,cz))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(ax,ay,az))
            w_POSCAR.write('%23.16f %21.16f %21.16f\n'%(bx,by,bz))
            for l_wPOSCAR in range(5,7+wn+is_Sd):
                w_POSCAR.write(allline[l_wPOSCAR]+'\n')
            for l_wPOSCAR in range(n_atom):
                w_POSCAR.write('%20.16f %19.16f %19.16f\n'%(newpos[l_wPOSCAR][2],newpos[l_wPOSCAR][0],newpos[l_wPOSCAR][1]))
            w_POSCAR.close()
            os.system('mv BPOSCAR oldPOSCAR')
            os.system('mv POSCARw BPOSCAR')


w_KPOINTS = open('KPOINTSw','w')
w_KPOINTS.write('KPATH\n'+str(knum)+'\nLine-Mode\nReciprocal\n')
len_kpath=len(kpath)
for i in range(len_kpath):
    w_KPOINTS.write('%12.8f%12.8f%12.8f    ! %-s\n'%(kpath[i][0][0],kpath[i][0][1],kpath[i][0][2],kmark[i][0]))
    w_KPOINTS.write('%12.8f%12.8f%12.8f    ! %-s\n\n'%(kpath[i][1][0],kpath[i][1][1],kpath[i][1][2],kmark[i][1]))
w_KPOINTS.close()

w_K = open('all-KPOINTS','w')
len_kpath=len(kpath)
for i in range(len_kpath):
    for j in range(knum):
        w_K.write('%12.8f%12.8f%12.8f 0.00\n'%(j*(kpath[i][1][0]-kpath[i][0][0])/(knum-1)+kpath[i][0][0],j*(kpath[i][1][1]-kpath[i][0][1])/(knum-1)+kpath[i][0][1],j*(kpath[i][1][2]-kpath[i][0][2])/(knum-1)+kpath[i][0][2]))
    # w_K.write('%12.8f%12.8f%12.8f    ! %-s\n\n'%(kpath[i][1][0],kpath[i][1][1],kpath[i][1][2],kmark[i][1]))
w_K.close()



f = open('BPOSCAR','r+')
l = f.read()
allline = l.split("\n")
f.close()
w_POSCAR = open('POSCARw','w')
for l_wPOSCAR in range(5):
    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
for l_wPOSCAR in range(5,7+wn+is_Sd):
    w_POSCAR.write(allline[l_wPOSCAR]+'\n')
for l_wPOSCAR in range(n_atom):
    w_POSCAR.write(allline[7+wn+is_Sd+l_wPOSCAR]+'\n')
w_POSCAR.close()




f = open('POSCARw','r+')
l = f.read()
allline = l.split("\n")
name = allline[0].split()[0]
ax,ay,az = [ float(x) for x in allline[2].split() ]
bx,by,bz = [ float(x) for x in allline[3].split() ]
cx,cy,cz = [ float(x) for x in allline[4].split() ]
l_n_atom = [ int(x) for x in allline[5+wn].split() ]
# l_atom = allline[0].split()
n_atom = sum(l_n_atom)
pos = []
is_Sd = 0 if len(allline[6+wn].split())==1 else 1
for l_pos in range(7+wn+is_Sd, 7+wn+n_atom+is_Sd):
    pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split() ]
    pos.append([pos_a,pos_b,pos_c])
f.close()
# n_element = len(l_atom)
# formula = ''.join([l_atom[i_element]+str(l_n_atom[i_element]) for i_element in range(n_element)])
#deal with structure for start
ra = [ax,ay,az]
rb = [bx,by,bz]
rc = [cx,cy,cz]

d=0.00001
mian=[]
# ra=[4,0,0]
# rb=[-1,3**0.5,0]
# rc=[0,0,1]
# ra = [-2.9212500000000001,4.2082499999999996,5.1470000000000002]
# rb = [2.9212500000000001,-4.2082499999999996,5.1470000000000002]
# rc = [2.9212500000000001,4.2082499999999996,-5.1470000000000002]
# ra=[ 3.2115200472080403,0,0]
# rb=[0,8.5240904989406658,0]
# rc=[0,0,4.5108099844478842]
ka=axn(cross(rb,rc),2*math.pi/(dot(ra,cross(rb,rc))))
kb=axn(cross(rc,ra),2*math.pi/(dot(rb,cross(rc,ra))))
kc=axn(cross(ra,rb),2*math.pi/(dot(rc,cross(ra,rb))))
b1=mat(ka)
b2=mat(kb)
b3=mat(kc)
# b1=mat([-1,1,1])
# b2=mat([1,-1,1])
# b3=mat([1,1,-1])
for i in np.arange(-1,2,1):
    for j in np.arange(-1,2,1):
        for k in np.arange(-1,2,1):
            if i==j==k==0:
                continue
            n=i*b1+j*b2+k*b3
            P=0.5*n[0,0]**2+0.5*n[0,1]**2+0.5*n[0,2]**2
            mian.append([n[0,0],n[0,1],n[0,2],P,1])

print(mian)
print(len(mian))
mian0=[]
for i in range(len(mian)-1):
    for j in range(i+1):
        mian1=mian[i+1]
        mian2=mian[j]
        if mian2[4]==0:
            continue
        if (0.5*mian1[0]*mian2[0]+0.5*mian1[1]*mian2[1]+0.5*mian1[2]*mian2[2]-mian2[3])\
            *(-mian2[3])<=d:
            mian[i+1][4]=0
        elif (0.5*mian2[0]*mian1[0]+0.5*mian2[1]*mian1[1]+0.5*mian2[2]*mian1[2]-mian1[3])\
            *(-mian1[3])<=d:
            mian[j][4]=0
for i in range(len(mian)):
    if mian[i][4]==1:
        mian0.append(mian[i])
print(mian0)
print(len(mian0))

xian=[]
for i in range(len(mian0)-1):
    for j in range(i+1):
        mian1=mian0[i+1]
        mian2=mian0[j]
        dian=[]
        for k in range(len(mian0)):
            mian3=mian0[k]
            C=mat([[mian1[0],mian1[1],mian1[2]],[mian2[0],mian2[1],mian2[2]],[mian3[0],mian3[1],mian3[2]]])
            P=mat([[mian1[3],mian2[3],mian3[3]]]).T
            if np.linalg.det(C)==0:
                continue
            X=C.I*P
            isin=0
            for l in range(len(mian0)):
                mian4=mian0[l]
                if (X[0,0]*mian4[0]+X[1,0]*mian4[1]+X[2,0]*mian4[2]-mian4[3])*(-mian4[3])<-d:
                    isin+=1
            if isin==0 and [X[0,0],X[1,0],X[2,0]] not in dian:
                dian.append([X[0,0],X[1,0],X[2,0]])
        if len(dian)>1 :
            xian.append(dian)
print(xian)
print(len(xian))
#fig = plt.figure(2)
# ax = Axes3D(fig)
#ax=fig.add_subplot(111,projection='3d')

#ax.plot([0,b1[0,0]],[0,b1[0,1]],[0,b1[0,2]],'b-')
#ax.plot([0,b2[0,0]],[0,b2[0,1]],[0,b2[0,2]],'b-')
#ax.plot([0,b3[0,0]],[0,b3[0,1]],[0,b3[0,2]],'b-')
#ax.text(b1[0,0],b1[0,1],b1[0,2],'ka',fontsize=15)
#ax.text(b2[0,0],b2[0,1],b2[0,2],'kb',fontsize=15)
#ax.text(b3[0,0],b3[0,1],b3[0,2],'kc',fontsize=15)

#for i in range(len(xian)):
#    print(i)
#    x=[xian[i][0][0],xian[i][1][0]]
#    y=[xian[i][0][1],xian[i][1][1]]
#    z=[xian[i][0][2],xian[i][1][2]]
#    print(x)
#    print(y)
#    print(z)
#    ax.plot(x,y,z,'k-')

#if os.path.exists('KPOINTSw'):
#    f = open('KPOINTSw','r+')
#    l = f.read().split('\n')
#    knum = int(l[1].split()[0])  #num of k point ver path
#    f.close()
#    lenl=len(l)
#    ib=0
#    signk=[]
#    while True:
#        if ib*3+6>lenl:
#            break
#        a=[]
#        signk.append(a)
#        signk[ib].append(l[ib*3+4].split()[4])
#        signk[ib].append(l[ib*3+5].split()[4])
#        ib=ib+1
#    nb=ib
#    f= open('KPOINTSw','r')
#    l=f.read()
#    f.close()
#    allline = l.split('\n')
#    for i in range(nb):
#        line1=allline[i*3+4].split()
#        x1=float(line1[0])
#        y1=float(line1[1])
#        z1=float(line1[2])
#        mark1=line1[4]
#        xian1=x1*b1+y1*b2+z1*b3
#
#        line2=allline[i*3+5].split()
#        x2=float(line2[0])
#        y2=float(line2[1])
#        z2=float(line2[2])
#        mark2=line2[4]
#        xian2=x2*b1+y2*b2+z2*b3

#        x=[xian1[0,0],xian2[0,0]]
#        y=[xian1[0,1],xian2[0,1]]
#        z=[xian1[0,2],xian2[0,2]]
#        ax.plot(x,y,z,'r')
#        ax.text(xian1[0,0],xian1[0,1],xian1[0,2],mark1,fontsize=15)
#        ax.text(xian2[0,0],xian2[0,1],xian2[0,2],mark2,fontsize=15)




#ax.set_xlabel('kx')
#ax.set_ylabel('ky')
#ax.set_zlabel('kz')

#ax.axis('scaled')
#ax.grid(False)
#plt.show()
#plt.savefig("brillouin.pdf")

#f=open('brillouin','w')
#for i in range(len(xian)):
#    f.write('%21.16f %21.16f %21.16f \n'%(xian[i][0][0],xian[i][0][1],xian[i][0][2]))
#    f.write('%21.16f %21.16f %21.16f \n'%(xian[i][1][0],xian[i][1][1],xian[i][1][2]))
    # f.write('\n')
#f.close()

print('KPATH writen in KPOINTSw, structure for band-calculation writen in POSCARw, brillouin plot in brillouin.pdf ,brillouin edges writen in brillouin, all the KPOINTS of band writen in all-KPOINTS')
