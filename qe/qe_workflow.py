from asyncio.log import logger
import os
import logging

from qe_inputpara import qe_inputpara
from qe_writesubmit import qe_writesubmit
from qe_writeinput import qe_writeinput
from qe_submitjob import qe_submitjob 

logger = logging.getLogger("qe_workflow")

class qe_workflow:
    def __init__(
        self,
        input_file_path: str, 
        work_path: str,
        **kwargs: dict,
    ):
        self.input_file_path     = input_file_path
        self.work_path           = work_path
        self.pressure            = None
        self.work_underpressure  = None
        self.kpoints_dense       = [None, None, None]
        self.kpoints_sparse      = [None, None, None] 
        self.qpoints             = [None, None, None]
        self.run_mode            = None
        self.submit_job_system   = "slurm"
        self.dyn0_flag           = False 
        self.q_non_irreducible_amount = None
        self.inserted_points_num = None

        if kwargs:
            for key, value in kwargs.items():
                if key == "pressure":
                    self.pressure = value
                    self.work_underpressure = os.path.join(self.work_path, str(self.pressure))
                if key == "kpoints_dense":
                    self.kpoints_dense  = value
                    self.kpoints_sparse = [kp/2 for kp in self.kpoints_dense]
                    self.qpoints        = [kp/4 for kp in self.kpoints_dense]
                if key == "kpoints_sparse":
                    self.kpoints_sparse = value
                    self.kpoints_dense  = [kp*2 for kp in self.kpoints_sparse]
                    self.qpoints        = [kp/2 for kp in self.kpoints_sparse]
                if key == "qpoints":
                    self.qpoints        = value 
                    self.kpoints_dense  = [kp*4 for kp in self.qpoints]
                    self.kpoints_sparse = [kp*2 for kp in self.qpoints]
                if key == "run_mode":
                    self.run_mode = value
                if key == "submit_job_system":
                    self.submit_job_system = value
                if key == "dyn0_flag":
                    self.dyn0_flag = value 
                if key== "q_non_irreducible_amount":
                    self.q_non_irreducible_amount = value
                if key== "inserted_points_num":
                    self.inserted_points_num = value

        if self.work_underpressure is None:
            self.work_underpressure = self.work_path 
        ############################ prepare pp directory #########################
        logger.info(f"create pp dir in {self.work_underpressure}")
        if self.input_file_path is None:
            logger.error("please specify the inputfile *.vasp or POSCAR")
            raise ValueError ("please specify the inputfile *.vasp or POSCAR")
        self.workpath_pppath = os.path.abspath(os.path.join(self.work_underpressure, "pp"))
        if not os.path.exists(self.workpath_pppath):
            os.makedirs(self.workpath_pppath)    
        ############################# done pp directory ##########################

        self.qe_inputpara = qe_inputpara(
            input_file_path=self.input_file_path,
            work_underpressure=self.work_underpressure,
            pressure=self.pressure,
            workpath_pppath=self.workpath_pppath,
            kpoints_dense=self.kpoints_dense,
            kpoints_sparse=self.kpoints_sparse,
            qpoints=self.qpoints,
            inserted_points_num=self.inserted_points_num,
            run_mode=self.run_mode,
        )

        self.qe_writeinput = qe_writeinput(
            qe_input_object=self.qe_inputpara, 
            run_mode=self.run_mode,
            q_non_irreducible_amount=self.q_non_irreducible_amount 
        )
        
        self.qe_writesubmit = qe_writesubmit(
            qe_input_object=self.qe_inputpara, 
            run_mode=self.run_mode,
            submit_job_system=self.submit_job_system,
            q_non_irreducible_amount=self.q_non_irreducible_amount 
        )

        self.qe_submitjob = qe_submitjob(
            qe_input_object=self.qe_inputpara, 
            run_mode=self.run_mode,
            submit_job_system=self.submit_job_system,
            dyn0_flag=self.dyn0_flag,
            q_non_irreducible_amount=self.q_non_irreducible_amount 
        )
        