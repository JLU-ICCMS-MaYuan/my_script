import os

from vasp.vasp_inputpara import vasp_inputpara 

class vasp_writeincar:
    
    def __init__(self, work_path, **kwargs) -> None:

        self.work_path = work_path
        for key, value in kwargs.items():
            setattr(self, key, value)

    def writeinput(self, mode=None, incar_path=None):
        if mode == None:
            mode = self.mode
        if incar_path == None:
            incar_path = self.work_path
        if mode == 'rvf' or mode == 'rv1':
            incar_dirpath = self.opt_fine_incar(incar_path)
        if mode == 'rv3':
            incar_dirpath1 = self.opt_incar1(incar_path)
            incar_dirpath2 = self.opt_incar2(incar_path)
            incar_dirpath  = self.opt_incar3(incar_path)
        if mode == 'rv4':
            incar_dirpath1 = self.opt_incar1(incar_path)
            incar_dirpath2 = self.opt_incar2(incar_path)
            incar_dirpath3 = self.opt_incar3(incar_path)
            incar_dirpath  = self.opt_incar4(incar_path)
        if mode == 'disp':
            incar_dirpath = self.disp_incar(incar_path)
        if mode == 'dfpt':
            incar_dirpath = self.dfpt_incar(incar_path)
        if mode == 'scf':
            incar_dirpath = self.scf_incar(incar_path)
        if mode == 'eband':
            incar_dirpath = self.eband_incar(incar_path)
        if mode == 'eledos':
            incar_dirpath = self.eledos_incar(incar_path)
        if mode == 'cohp':
            incar_dirpath = self.cohp_incar(incar_path)

        if int(self.ispin)== 2:
            self.append_magnet(incar_dirpath)
        if self.ldau  == ".TRUE.":
            self.append_u(incar_dirpath)

    @classmethod
    def init_from_relaxinput(cls, other_class: vasp_inputpara):
        self = cls(
            work_path=other_class.work_path,
            encut=other_class.encut,
            kspacing=other_class.kspacing,
            ismear=other_class.ismear,
            sigma=other_class.sigma,
            ediff=other_class.ediff,
            ediffg=other_class.ediffg,
            ibrion=other_class.ibrion,
            isif=other_class.isif,
            potim=other_class.potim,
            nelm=other_class.nelm,
            ncore=other_class.ncore,
            lreal=other_class.lreal,
            symprec=other_class.symprec,
            
            press=other_class.press,
            mode=other_class.mode,

            # consider magnetism
            isym=other_class.isym,
            ispin=other_class.ispin,
            magmom=other_class.magmom,
            lorbit=other_class.lorbit,
            lasph=other_class.lasph,
            gga=other_class.gga,

            # DFT+U+J
            ldau=other_class.ldau,
            ldautype=other_class.ldautype,
            ldaul=other_class.ldaul,
            ldauu=other_class.ldauu,
            ldauj=other_class.ldauj,
            lmaxmix=other_class.lmaxmix,
        )
        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: vasp_inputpara):
        self = cls(
            work_path=other_class.work_path,
            encut=other_class.encut,
            ismear=other_class.ismear,
            sigma=other_class.sigma,
            ediff=other_class.ediff,
            ediffg=other_class.ediffg,
            potim=other_class.potim,
            nelm=other_class.nelm,
            ncore=other_class.ncore,
            lreal=other_class.lreal,
            mode=other_class.mode,
            kspacing=other_class.kspacing,
            symprec=other_class.symprec,

            # consider magnetism
            isym=other_class.isym,
            ispin=other_class.ispin,
            magmom=other_class.magmom,
            lorbit=other_class.lorbit,
            lasph=other_class.lasph,
            gga=other_class.gga,

            # DFT+U+J
            ldau=other_class.ldau,
            ldautype=other_class.ldautype,
            ldaul=other_class.ldaul,
            ldauu=other_class.ldauu,
            ldauj=other_class.ldauj,
            lmaxmix=other_class.lmaxmix,
        )
        return self

    @classmethod
    def init_from_eletron(cls, other_class: vasp_inputpara):
        self = cls(
            work_path=other_class.work_path,
            encut=other_class.encut,
            ismear=other_class.ismear,
            sigma=other_class.sigma,
            ediff=other_class.ediff,
            ediffg=other_class.ediffg,
            nelm=other_class.nelm,
            ncore=other_class.ncore,
            lreal=other_class.lreal,
            mode=other_class.mode,
            kspacing=other_class.kspacing,
            npar=other_class.npar,
            nbands=other_class.nbands,

            # consider magnetism
            isym=other_class.isym,
            ispin=other_class.ispin,
            magmom=other_class.magmom,
            lorbit=other_class.lorbit,
            lasph=other_class.lasph,
            gga=other_class.gga,

            # DFT+U+J
            ldau=other_class.ldau,
            ldautype=other_class.ldautype,
            ldaul=other_class.ldaul,
            ldauu=other_class.ldauu,
            ldauj=other_class.ldauj,
            lmaxmix=other_class.lmaxmix,

            # Eledos
            nedos=other_class.nedos
        )
        return self

    def opt_incar1(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_1")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")
            incar.write("ISYM     = {}   \n".format(str(self.isym))) 
            incar.write("ENCUT    = 350  \n")
            incar.write("PREC     = LOW  \n") 
            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.8  \n")            
            incar.write("ISMEAR   = 0    \n")   
            incar.write("SIGMA    = 0.2  \n")   
            incar.write("NELM     = 200   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-3 \n")
 
            incar.write("NSW      = 200  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 2    \n")
            incar.write("POTIM    = 0.3  \n")
            incar.write("LWAVE  = .FALSE.\n")                 
            incar.write("LCHARG = .FALSE.\n")   
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))
        return incar_filepath

    def opt_incar2(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_2")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")  
            incar.write("ISYM     = {}   \n".format(str(self.isym))) 
            incar.write("ENCUT    = 400  \n")        
            incar.write("PREC     = Normal\n") 
            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.5  \n")
            incar.write("ISMEAR   = 0    \n")   
            incar.write("SIGMA    = 0.2  \n")   
            incar.write("NELM     = 200   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-4 \n")
            incar.write("EDIFFG   = -0.2 \n")  
            incar.write("NSW      = 400  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 4    \n")   
            incar.write("POTIM    = 0.1  \n")
            incar.write("LWAVE  = .FALSE.\n")                 
            incar.write("LCHARG = .FALSE.\n")   
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)) )  
        return incar_filepath

    def opt_incar3(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_3")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")
            incar.write("ISYM     = {}   \n".format(str(self.isym))) 
            incar.write("ENCUT    = 500  \n")        
            incar.write("PREC     = A    \n")
            incar.write("SYMPREC  = {}   \n".format(str(self.symprec)))

            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.30 \n")            
            incar.write("ISMEAR   = 0    \n")   
            incar.write("SIGMA    = 0.05 \n")   
            incar.write("NELM     = 90   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-5 \n")
            incar.write("EDIFFG   = -0.1\n")
            incar.write("NSW      = 500  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 3    \n")        
            incar.write("POTIM    = 0.05 \n")
            incar.write("LWAVE    = .FALSE.\n")                 
            incar.write("LCHARG   = .FALSE.\n")   
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10))) 
        return incar_filepath   

    def opt_incar4(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_4")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")
            incar.write("ISYM     = {}   \n".format(str(self.isym))) 
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
            incar.write("SYMPREC  = {}   \n".format(str(self.symprec)))

            incar.write("NCORE    = {}   \n".format(str(self.ncore)))         
            incar.write("KSPACING = {}   \n".format(str(self.kspacing))) 
            incar.write("ISMEAR   = {}   \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}   \n".format(str(self.sigma)))   
            incar.write("NELM     = {}   \n".format(str(self.nelm)))   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = {}   \n".format(str(self.ediff)))
            incar.write("EDIFFG   = {}   \n".format(str(self.ediffg)))
            incar.write("NSW      = 500  \n")   
            incar.write("IBRION   = {}   \n".format(str(self.ibrion)))   
            incar.write("ISIF     = {}   \n".format(str(self.isif)))    
            incar.write("POTIM    = {}   \n".format(str(self.potim)))
            incar.write("LWAVE  = .FALSE.\n")                 
            incar.write("LCHARG = .FALSE.\n")   
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))  
        return incar_filepath 

    def opt_fine_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")  
            incar.write("ISYM     = {}   \n".format(str(self.isym))) 
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
            incar.write("SYMPREC  = {}   \n".format(str(self.symprec)))

            incar.write("NCORE    = {}   \n".format(str(self.ncore)))         
            incar.write("KSPACING = {}   \n".format(str(self.kspacing))) 
            incar.write("ISMEAR   = {}   \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}   \n".format(str(self.sigma)))   
            incar.write("NELM     = {}   \n".format(str(self.nelm)))   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = {}   \n".format(str(self.ediff)))
            incar.write("EDIFFG   = {}   \n".format(str(self.ediffg)))
            incar.write("NSW      = 500  \n")   
            incar.write("IBRION   = {}   \n".format(str(self.ibrion)))   
            incar.write("ISIF     = {}   \n".format(str(self.isif)))    
            incar.write("POTIM    = {}   \n".format(str(self.potim)))
            incar.write("LWAVE  = .FALSE.\n")                 
            incar.write("LCHARG = .FALSE.\n")           

            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))  
        return incar_filepath
            
    def disp_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
            incar.write("SYMPREC  = {}   \n".format(str(self.symprec)))
            incar.write("ISMEAR   = {}   \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}   \n".format(str(self.sigma)))   
            incar.write("NELM     = {}   \n".format(str(self.nelm)))   
            incar.write("EDIFF    = {}   \n".format(self.ediff))
            incar.write("EDIFFG   = {}   \n".format(self.ediffg))
            incar.write("IBRION   = -1   \n")   
            incar.write("IALGO    = 38   \n")

            incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL    = {}    \n".format(str(self.lreal)))
            incar.write("LWAVE    =.FALSE.\n")
            incar.write("LCHARG   =.FALSE.\n")
            incar.write("ADDGRID  = .TRUE.\n")
            if self.kspacing is not None:
                incar.write("KSPACING = {}   \n".format(str(self.kspacing)))
        return incar_filepath

    def dfpt_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0      \n")   
            incar.write("ICHARG   = 2      \n")   
            incar.write("ENCUT    = {}     \n".format(str(self.encut)))        
            incar.write("PREC     = A      \n")
            incar.write("SYMPREC  = {}   \n".format(str(self.symprec)))
            incar.write("ISMEAR   = {}     \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}     \n".format(str(self.sigma)))   
            incar.write("NELM     = {}     \n".format(str(self.nelm)))   
            incar.write("NELMIN   = 6      \n")   
            incar.write("EDIFF    = {}     \n".format(self.ediff))
            incar.write("EDIFFG   = {}     \n".format(self.ediffg))
            incar.write("IBRION   = 8      \n")   
            incar.write("IALGO    = 38     \n")
            incar.write("POTIM    = 0.01   \n") 

            incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL    = {}    \n".format(str(self.lreal)))
            incar.write("LWAVE    = .FALSE.\n")  
            incar.write("LCHARG   = .FALSE.\n") 
            incar.write("ADDGRID  = .TRUE. \n")
            if self.kspacing is not None:
                incar.write("KSPACING = {}   \n".format(str(self.kspacing)))
        return incar_filepath

    def scf_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0      \n")   
            incar.write("ICHARG   = 2      \n")  
            incar.write("ISYM     = {}     \n".format(str(self.isym))) 
            incar.write("ENCUT    = {}     \n".format(self.encut))       
            incar.write("PREC     = Accurate\n")
            incar.write("ISMEAR   = {}     \n".format(self.ismear))
            incar.write("SIGMA    = {}     \n".format(self.sigma)) 
            incar.write("NELM     = {}     \n".format(self.nelm))
            incar.write("NELMIN   = 2      \n")   
            incar.write("EDIFF    = {}     \n".format(self.ediff))
            incar.write("EDIFFG   = {}     \n".format(self.ediffg))
            incar.write("IBRION   = -1     \n")   
            incar.write("NSW      = 0      \n")
            incar.write("KSPACING = {}     \n".format(self.kspacing))
            incar.write("KGAMMA   = .TRUE. \n") # Determines whether the k points (specified by the KSPACING tag ) include (KGAMMA=.TRUE.) the Γ\Gamma  point.
            incar.write("VOSKOWN  = 1      \n")
            incar.write("NBLOCK   = 1      \n")
            incar.write("NWRITE   = 1      \n")
            incar.write("ALGO     = Normal \n")
            incar.write("ISPIN    = 1      \n")
            incar.write("INIWAV   = 1      \n")
            incar.write("NBANDS   = {}     \n".format(self.nbands))
            #incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL  = {}       \n".format(self.lreal))
            incar.write("LWAVE    = .TRUE.  \n")  
            incar.write("ADDGRID  = .TRUE.  \n")
            incar.write("#RWIGS   = 1.54 0.82\n")
            incar.write("LHYPERFINE = .FALSE.\n")
            incar.write("NPAR   = {}         \n".format(self.npar))          
            incar.write("\n")
            incar.write("LCHARG  =.TRUE.    # For Bader  \n") # 能带计算需要将其打开  确保这个是TRUE            
            incar.write("LAECHG  =.TRUE.    # For Bader  \n")   
            incar.write("LELF    =.TRUE.    # For ELF    \n")
        return incar_filepath

    def eband_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART = 0            \n")
            incar.write("ICHARG = 11           \n")
            incar.write("ISYM     = {}         \n".format(str(self.isym))) 
            incar.write("ENCUT  = {}           \n".format(self.encut))       
            incar.write("PREC   = Accurate     \n") 
            incar.write("NELM   = {}           \n".format(self.nelm))
            incar.write("ISMEAR = {}           \n".format(self.ismear))
            incar.write("SIGMA  = {}           \n".format(self.sigma)) 
            incar.write("NELMIN = 2            \n")
            incar.write("EDIFF  = {}           \n".format(self.ediff))
            incar.write("IBRION = -1           \n")             
            incar.write("NSW    = 0            \n")                         
            incar.write("VOSKOWN= 1            \n")       
            incar.write("NWRITE = 3            \n")            
            incar.write("ALGO   = Normal       \n")                                  
            incar.write("ISPIN  = 1            \n")           
            incar.write("INIWAV = 1            \n")            
            incar.write("LREAL  = {}           \n".format(self.lreal))
            incar.write("NBANDS = {}           \n".format(self.nbands))
            incar.write("LWAVE  = .TRUE.      \n")                 
            incar.write("LCHARG = .FALSE.      \n")           
            incar.write("ADDGRID= .TRUE.       \n")   
            incar.write("#RWIGS = 1.54 0.82    \n")     
            incar.write("LHYPERFINE = .FALSE.  \n")                      
            incar.write("NPAR   = {}           \n".format(self.npar))          
            incar.write("LORBIT = 11           \n") # 算投影能带有用
        return incar_filepath

    def eledos_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART  = 1            \n") # if a WAVECAR file exists
            incar.write("ICHARG  = 11           \n") # 从CHGCAR读取给定电荷密度的特征值(用于带结构图)或状态密度(DOS)。自洽CHGCAR文件必须事先通过一个跨越整个布里渊区的k点网格进行完全自洽计算来确定。
            incar.write("ISYM    = {}           \n".format(str(self.isym))) 
            incar.write("ENCUT   = {}           \n".format(self.encut))           
            incar.write("PREC    = Accurate     \n") 
            incar.write("ISMEAR  = -5      # For DOS\n")
            incar.write("SIGMA   = {}           \n".format(self.sigma))
            incar.write("NELM    = 100          \n")                               
            incar.write("NELMIN  = 2            \n")
            incar.write("EDIFF   = {}           \n".format(self.ediff))
            incar.write("IBRION  = -1           \n")             
            incar.write("NSW     = 0            \n")                         
            incar.write("VOSKOWN = 1            \n")       
            incar.write("NWRITE  = 3            \n")            
            incar.write("ALGO    = Normal       \n")                                  
            incar.write("ISPIN   = 1            \n")           
            incar.write("INIWAV  = 1            \n")            
            incar.write("LREAL   = {}           \n".format(self.lreal))
            incar.write("NBANDS  = {}           \n".format(self.nbands))
            incar.write("LWAVE   = .FALSE.      \n")                 
            incar.write("ADDGRID = .TRUE.       \n")   
            incar.write("#RWIGS  = 1.54 0.82    \n")     
            incar.write("LHYPERFINE = .FALSE.   \n")                      
            incar.write("NPAR    = {}           \n".format(self.npar))          
            incar.write("NEDOS   = {}           \n".format(self.nedos)) # NEDOS指定DOS被评估的网格点的数量
            incar.write("LORBIT  = 11           \n") # 输出分波态密度信息
            incar.write("#EMIN   = -10          \n") # 此为DOS图的能量范围，根据能带的能量范围来决定min和max是多少。
            incar.write("#EMAX   =  10          \n") 
            incar.write("LCHARG  = .TRUE.       \n")                  
        return incar_filepath

    def cohp_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART  = 0            \n") # if a WAVECAR file exists
            incar.write("ICHARG  = 2            \n") # 从CHGCAR读取给定电荷密度的特征值(用于带结构图)或状态密度(DOS)。自洽CHGCAR文件必须事先通过一个跨越整个布里渊区的k点网格进行完全自洽计算来确定。
            incar.write("ISYM    =-1            \n")
            incar.write("ENCUT   = {}           \n".format(self.encut))          
            incar.write("PREC    = Accurate     \n") 
            incar.write("ISMEAR  = -5           \n")
            incar.write("NELM    = {}           \n".format(self.nelm))                               
            incar.write("NELMIN  = 2            \n")
            incar.write("EDIFF   = {}           \n".format(self.ediff))
            incar.write("IBRION  = -1           \n")
            incar.write("ISIF    = 0            \n")
            incar.write("NSW     = 0            \n")
            incar.write("ISPIN   = 1            \n")
            incar.write("LREAL   = .FALSE.      \n")
            incar.write("NBANDS  = {}           \n".format(self.nbands))
            incar.write("NPAR    = {}           \n".format(self.npar))          
            incar.write("NEDOS   = {}           \n".format(self.nedos)) # NEDOS指定DOS被评估的网格点的数量
            incar.write("LORBIT  = 12           \n") # 输出分波态密度信息
            incar.write("#EMIN   = -10          \n") # 此为DOS图的能量范围，根据能带的能量范围来决定min和max是多少。
            incar.write("#EMAX   =  10          \n") 
            incar.write("LREAL   = {}           \n".format(self.lreal))
            incar.write("LWAVE   = .TRUE.       \n")  
            incar.write("ADDGRID = .TRUE.       \n")
        return incar_filepath

    def append_magnet(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_dirpath, "a") as incar:
            incar.write("\n")
            incar.write("# add parameters for magnetism\n")
            incar.write("ISYM     = {}\n".format(self.isym))
            incar.write("ISPIN    = {}\n".format(self.ispin))
            incar.write("MAGMOM   = {}\n".format(self.magmom))
            incar.write("LORBIT   = {}\n".format(self.lorbit))
            incar.write("LASPH    = {}\n".format(self.lasph))
            incar.write("GGA      = {}\n".format(self.gga))
              
    def append_u(self, incar_dirpath):
        with open(incar_dirpath, "a") as incar:
            incar.write("\n")
            incar.write("# add parameters for U value\n")
            incar.write("LDAU     = {}\n".format(self.ldau))
            incar.write("LDAUTYPE = {}\n".format(self.ldautype))
            incar.write("LDAUL    = {}\n".format(self.ldaul))
            incar.write("LDAUU    = {}\n".format(self.ldauu))
            incar.write("LDAUJ    = {}\n".format(self.ldauj)) 
            incar.write("LMAXMIX  = {}\n".format(self.lmaxmix))
