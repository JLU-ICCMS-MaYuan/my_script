from asyncio.log import logger
import os
import math

from pymatgen.io.vasp import Kpoints

from vasp_inputpara import vasp_inputpara


class vasp_writekpoints:
    def __init__(
        self, 
        input_file_path: str,
        work_underpressure: str,
        **kwargs, 
        ):
        self.input_file_path = input_file_path
        self.work_underpressure = work_underpressure
        self.output_kpoints = os.path.join(self.work_underpressure, "KPOINTS") 
        for key, value in kwargs.items():
            setattr(self, key, value)        

        if self.kpoints != [None, None, None]:
            logger.info(f"you set the kpoints number is {self.kpoints}")
            self.create_kpoints(self.output_kpoints)
        else:
            logger.info(f"you set the kspacing is {self.kspacing}")
            self.create_mesh(self.kspacing, self.input_file_path, self.output_kpoints)

    @classmethod
    def init_from_inputpara(cls, other_class: vasp_inputpara):
        self = cls(
            input_file_path=other_class.input_file_path,
            work_underpressure=other_class.work_underpressure,
            sposcar_struct_type=other_class.sposcar_struct_type,
            kpoints=other_class.kpoints,
            kspacing=other_class.kspacing,
        )
        return self

    def create_kpoints(self, output_kpoints):
        kpoints = Kpoints.automatic_density_by_lengths(
            self.sposcar_struct_type, 
            length_densities=self.kpoints,
            force_gamma=True)
        kpoints.write_file(output_kpoints)
        print(kpoints)

    def create_mesh(self, kresolution, input_file_path, output_kpoints):
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
        n_k1 = math.ceil(b1/interval_k)
        n_k2 = math.ceil(b2/interval_k)
        n_k3 = math.ceil(b3/interval_k)
        newpos = []
        for i_natom in range(n_atom):
            newpos.append([pos[i_natom][0],pos[i_natom][1],pos[i_natom][2]])

        #writw KPOINTS
        w_KPOINTS = open(output_kpoints, 'w')

        w_KPOINTS.write('Automatic mesh\n0\nGamma\n  %-3s %-3s %-3s\n  0.  0.  0.\n' %(int(n_k1),int(n_k2),int(n_k3)))
        w_KPOINTS.close()

        print('Kmesh is writen in KPOINTS-mesh')
