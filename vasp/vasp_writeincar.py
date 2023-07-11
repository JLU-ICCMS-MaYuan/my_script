import os

from vasp.vasp_inputpara import vasp_inputpara 

class vasp_writeincar:
    
    def __init__(self, work_path, **kwargs) -> None:

        self.work_path = work_path
        for key, value in kwargs.items():
            setattr(self, key, value)

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
        )
        return self

    def opt_incar1(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_1")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = 350  \n")
            incar.write("PREC     = LOW  \n") 
            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.8  \n")            
            incar.write("ISMEAR   = 1    \n")   
            incar.write("SIGMA    = 0.2  \n")   
            incar.write("NELM     = 200   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-3 \n")
 
            incar.write("NSW      = 200  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 2    \n")
            incar.write("POTIM    = 0.3  \n")
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))   

    def opt_incar2(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_2")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = 400  \n")        
            incar.write("PREC     = Normal\n") 
            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.5  \n")
            incar.write("ISMEAR   = 1    \n")   
            incar.write("SIGMA    = 0.2  \n")   
            incar.write("NELM     = 200   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-4 \n")
            incar.write("EDIFFG   = -0.2 \n")  
            incar.write("NSW      = 400  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 4    \n")   
            incar.write("POTIM    = 0.1  \n")
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)) )   

    def opt_incar3(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_3")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = 500  \n")        
            incar.write("PREC     = A    \n")
            incar.write("NCORE    = 4    \n")         
            incar.write("KSPACING = 0.30 \n")            
            incar.write("ISMEAR   = 1    \n")   
            incar.write("SIGMA    = 0.05 \n")   
            incar.write("NELM     = 90   \n")   
            incar.write("NELMIN   = 6    \n")   
            incar.write("EDIFF    = 1e-7 \n")
            incar.write("EDIFFG   = -0.001\n")
            incar.write("NSW      = 500  \n")   
            incar.write("IBRION   = 2    \n")   
            incar.write("ISIF     = 3    \n")        
            incar.write("POTIM    = 0.05 \n")
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))    

    def opt_fine_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
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
            incar.write("PSTRESS  = {}   \n".format(str(float(self.press)*10)))  

            
    def disp_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
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


    def dfpt_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0      \n")   
            incar.write("ICHARG   = 2      \n")   
            incar.write("ENCUT    = {}     \n".format(str(self.encut)))        
            incar.write("PREC     = A      \n")
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

                
    def scf_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0      \n")   
            incar.write("ICHARG   = 2      \n")   
            incar.write("ENCUT    = {}     \n".format(str(self.encut)))        
            incar.write("PREC     = Accurate\n")
            incar.write("ISMEAR   = {}     \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}     \n".format(str(self.sigma)))   
            incar.write("NELM     = {}     \n".format(str(self.nelm)))   
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
            incar.write("#NBANDS   = 64     \n")
            #incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL    = .FALSE. \n")
            incar.write("LWAVE    = .FALSE.\n")  
            incar.write("LCHARG   = .TRUE. \n")  # 能带计算需要将其打开  确保这个是TRUE   
            incar.write("ADDGRID  = .FALSE.\n")
            incar.write("#RWIGS   = 1.54 0.82\n")
            incar.write("LHYPERFINE = .FALSE.\n")
            incar.write("NPAR     = 4      \n")

    def band_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART = 0           \n")
            incar.write("ICHARG = 11          \n")
            incar.write("ENCUT = 850          \n")             
            incar.write("PREC = Accurate      \n") 
            incar.write("NELM = 100           \n")                               
            incar.write("NELMIN = 2           \n")
            incar.write("EDIFF = 1.0e-08      \n")         
            incar.write("IBRION = -1          \n")             
            incar.write("NSW = 0              \n")                         
            incar.write("VOSKOWN = 1          \n")       
            incar.write("NWRITE = 3           \n")            
            incar.write("ALGO = Normal        \n")                                  
            incar.write("ISPIN = 1            \n")           
            incar.write("INIWAV = 1           \n")            
            incar.write("LREAL = .FALSE.      \n")
            incar.write("#NBANDS = 64         \n")              
            incar.write("LWAVE = .FALSE.      \n")                 
            incar.write("LCHARG = .FALSE.      \n")           
            incar.write("ADDGRID = .FALSE.    \n")   
            incar.write("#RWIGS = 1.54 0.82   \n")     
            incar.write("LHYPERFINE = .FALSE. \n")                      
            incar.write("NPAR = 4             \n")          
            incar.write("LORBIT = 11          \n") # 算投影能带有用


    def eledos_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART = 1           \n") # if a WAVECAR file exists
            incar.write("ICHARG = 11          \n") # 从CHGCAR读取给定电荷密度的特征值(用于带结构图)或状态密度(DOS)。自洽CHGCAR文件必须事先通过一个跨越整个布里渊区的k点网格进行完全自洽计算来确定。
            incar.write("ENCUT  = {}        \n".format(str(self.encut)))          
            incar.write("PREC = Accurate      \n") 
            incar.write("ISMEAR  = -5         \n")
            incar.write("SIGMA   = {}     \n".format(str(self.sigma)))   
            incar.write("NELM = 100           \n")                               
            incar.write("NELMIN = 2           \n")
            incar.write("EDIFF    = {}     \n".format(self.ediff))
            incar.write("IBRION = -1          \n")             
            incar.write("NSW = 0              \n")                         
            incar.write("VOSKOWN = 1          \n")       
            incar.write("NWRITE = 3           \n")            
            incar.write("ALGO = Normal        \n")                                  
            incar.write("ISPIN = 1            \n")           
            incar.write("INIWAV = 1           \n")            
            incar.write("LREAL = .FALSE.      \n")
            incar.write("#NBANDS = 64         \n")              
            incar.write("LWAVE = .FALSE.      \n")                 
            incar.write("LCHARG =.FALSE.      \n")                  
            incar.write("ADDGRID = .FALSE.    \n")   
            incar.write("#RWIGS = 1.54 0.82   \n")     
            incar.write("LHYPERFINE = .FALSE. \n")                      
            incar.write("NPAR = 4             \n") #这个应该是DOS的采点个数，弄高一点无所谓。
            incar.write("NEDOS = 1201         \n") # NEDOS指定DOS被评估的网格点的数量
            incar.write("LORBIT = 11          \n") # 输出分波态密度信息
            incar.write("#EMIN = -10           \n") # 此为DOS图的能量范围，根据能带的能量范围来决定min和max是多少。
            incar.write("#EMAX =  10           \n") 
