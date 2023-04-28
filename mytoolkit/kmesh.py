#!/usr/bin/env python
import math

def create_kmesh(kresolution, input_file_path, output_kpoints):
    """
    input parameter:
        kresolution: the meshing density of kpoints, equal to kspacing
        input_file_path: the path of poscar
        output_kpoints: output the file of KPOINTS
    """
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

    f = open(input_file_path, 'r+')
    l = f.read()
    allline = l.split("\n")
    name = allline[0].split()[0]
    ax,ay,az = [ float(x) for x in allline[2].split() ]
    bx,by,bz = [ float(x) for x in allline[3].split() ]
    cx,cy,cz = [ float(x) for x in allline[4].split() ]
    l_n_atom = [ int(x) for x in allline[6].split() ]
    l_atom = allline[5].split()
    n_atom = sum(l_n_atom)
    pos = []
    is_Sd = 0 if len(allline[7].split())==1 else 1
    for l_pos in range(8+is_Sd, 8+n_atom+is_Sd):
        pos_a,pos_b,pos_c = [ float(x) for x in allline[l_pos].split()[0:3] ]
        pos.append([pos_a,pos_b,pos_c])
    f.close()
    n_element = len(l_atom)
    formula = ''.join([l_atom[i_element]+str(l_n_atom[i_element]) for i_element in range(n_element)])
    #deal with structure for start
    interval_k = float(kresolution)
    a1 = [ax,ay,az]
    a2 = [bx,by,bz]
    a3 = [cx,cy,cz]
    a1mol=mol(a1)
    a2mol=mol(a2)
    a3mol=mol(a3)
    b1=axn(cross(a2,a3),2*math.pi/(dot(a1,cross(a2,a3))))
    b2=axn(cross(a3,a1),2*math.pi/(dot(a2,cross(a3,a1))))
    b3=axn(cross(a1,a2),2*math.pi/(dot(a3,cross(a1,a2))))
    b1mol=mol(b1)
    b2mol=mol(b2)
    b3mol=mol(b3)
    n_k1 = math.ceil(b1mol/interval_k)
    n_k2 = math.ceil(b2mol/interval_k)
    n_k3 = math.ceil(b3mol/interval_k)
    newpos = []
    for i_natom in range(n_atom):
        newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])

    #writw KPOINTS
    w_KPOINTS = open(output_kpoints, 'w')

    w_KPOINTS.write('Automatic mesh\n0\nGamma\n  %-3s %-3s %-3s\n  0.  0.  0.\n' %(int(n_k1),int(n_k2),int(n_k3)))
    w_KPOINTS.close()

    print("%-3s %-3s %-3s" % (int(n_k1),int(n_k2),int(n_k3)))
    print('Kmesh is writen in KPOINTS-mesh')


if __name__ == "__main__":
    kresolution = float(input("(注意: 这个脚本的作用时将KSPACING转化为对应的KPOINTS, \n所应用的公式就是vasp官网给出的: N_i = max(1, ceiling(|b_i|*2*pi/KSPACING)), \n这里需要你输入的k点密度不同于vaspkit中需要你输入的k点密度。)\n请输入k点密度: \n"))
    create_kmesh(kresolution, "POSCAR", "KPOINTS")