import os

from vasp_inputpara import vasp_inputpara 

class vasp_writeincar:
    
    def __init__(self, work_underpressure, **kwargs) -> None:

        self.work_underpressure = work_underpressure
        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.mode == 'rvf':
            self.opt_fine_incar(self.work_underpressure) 
        elif self.mode == 'rv3':
            self.opt_incar1(self.work_underpressure)
            self.opt_incar2(self.work_underpressure)
            self.opt_incar3(self.work_underpressure)
        elif self.mode == 'disp':
            self.disp_incar(self.work_underpressure)
        elif self.mode == 'dfpt':
            self.dfpt_incar(self.work_underpressure)

    @classmethod
    def init_from_relaxinput(cls, other_class: vasp_inputpara):
        self = cls(
            work_underpressure=other_class.work_underpressure,
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
            
            press=other_class.press,
            mode=other_class.mode,
        )
        return self

    @classmethod
    def init_from_phonoinput(cls, other_class: vasp_inputpara):
        self = cls(
            work_underpressure=other_class.work_underpressure,
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
        incar_filepath = os.path.join(incar_dirpath, "INCAR_fine")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
            
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
        incar_filepath = os.path.join(incar_dirpath, "INCAR_disp")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0    \n")   
            incar.write("ICHARG   = 2    \n")   
            incar.write("ENCUT    = {}   \n".format(str(self.encut)))        
            incar.write("PREC     = A    \n")
            incar.write("ISMEAR   = {}   \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}   \n".format(str(self.sigma)))   
            incar.write("NELM     = {}   \n".format(str(self.nelm)))   
            incar.write("EDIFF    = 1e-6 \n")
            incar.write("EDIFFG   =-0.01\n")
            incar.write("IBRION   = -1   \n")   
            incar.write("IALGO    = 38   \n")

            #incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL    = {}    \n".format(str(self.lreal)))
            incar.write("LWAVE    =.FALSE.\n")
            incar.write("LCHARG   =.FALSE.\n")
            incar.write("ADDGRID  = .TRUE.\n")

    def dfpt_incar(self, incar_dirpath):
        incar_filepath = os.path.join(incar_dirpath, "INCAR_dfpt")
        with open(incar_filepath, "w") as incar:
            incar.write("ISTART   = 0      \n")   
            incar.write("ICHARG   = 2      \n")   
            incar.write("ENCUT    = {}     \n".format(str(self.encut)))        
            incar.write("PREC     = A      \n")
            incar.write("ISMEAR   = {}     \n".format(str(self.ismear)))   
            incar.write("SIGMA    = {}     \n".format(str(self.sigma)))   
            incar.write("NELM     = {}     \n".format(str(self.nelm)))   
            incar.write("NELMIN   = 6      \n")   
            incar.write("EDIFF    = 1e-6   \n")
            incar.write("EDIFFG   = -0.01 \n")
            incar.write("IBRION   = 8      \n")   
            incar.write("IALGO    = 38     \n")
            incar.write("POTIM    = 0.01   \n") 

            #incar.write("NCORE    = {}    \n".format(str(self.ncore)))         
            incar.write("LREAL    = {}    \n".format(str(self.lreal)))
            incar.write("LWAVE    = .FALSE.\n")  
            incar.write("LCHARG   = .FALSE.\n") 
            incar.write("ADDGRID  = .TRUE. \n")









