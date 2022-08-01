#!/public/home/mayuan/miniconda3/envs/cage/bin/python3

"""
# check tasks number
j=0;x=1; for i in `squeue | awk '{print $1}'`; do  let j+=x; done; echo $j

relax:24 24 24
    qe_main.py -i ./test/POSCAR -w ./test/out -kd 24 24 24 -p 150 -mode relax -j pbs 

    qe_main.py -w ./ -scffit -kd 24 24 24 
    qe_main.py -w ./test/out/150.0/  -kd 16 16 16 -mode scffit -j pbs

    qe_main.py -w ./ -scf -ks 12 12 12 


    # no split method
    qe_main.py -w ./ -ph-nosplit -q 4 4 4


    # spilit method 1
    qe_main.py -w ./ -ph-nosplit -q 4 4 4 -dyn0
    qe_main.py -w ./ -ph-split-from-dyn0 -q 4 4 4
    qe_main.py -w ./ -m -q 4 4 4 

    # split method 2
    qe_main.py -w ./ -ph-split-startlastq -q 4 4 4 -q-no-irr-num 8

    # decide the system of Job submission system



"""
import os
import re
import time
import logging
from pathlib import Path
from argparse import ArgumentParser

from qe_workflow import qe_workflow

logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":


    logger.info("Start qe calculate")
    parser = ArgumentParser()
    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your POSCAR path, please notice your file format had better to be ***.vasp",
    )
    parser.add_argument(
        '-p',
        '--press',
        type=float,
        default=0.0,
        dest='press',
        help="please tell me your press which you were on",
    )
    parser.add_argument(
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory, I will put all file there",
    )
    parser.add_argument(
        '-kd',
        '-kpoints_dense',
        type=int,
        nargs="+",
        action='store',
        default=[16, 16, 16],
        dest='kpoints_dense',
        help="please tell me your kpoints_dense dense!",
    )
    parser.add_argument(
        '-ks',
        '-kpoints_sparse',
        type=int,
        nargs="+",
        action='store',
        default=[1, 1, 1],
        dest='kpoints_sparse',
        help="please tell me your kpoints_sparse dense!",
    )
    parser.add_argument(
        '-q',
        '-qpoints',
        type=int,
        nargs="+",
        action='store',
        default=[1, 1, 1],
        dest='qpoints',
        help="please tell me your qpoints dense!",
    )    
    parser.add_argument(
        '-mode',
        '--running-mode',
        type=str,
        action='store',
        default=None,
        dest='run_mode',
        help="please tell me your running mode!",
    )
    parser.add_argument(
        '-j',
        '-submit-job-system',
        type=str,
        default="slurm",
        dest='submit_job_system',
        help="please tell me your job submition system, eg: slurm, pbs",
    ) 
    args = parser.parse_args()
    input_file_path       = args.input_file_path
    press                 = args.press
    work_path             = args.work_path

    kpoints_dense     = args.kpoints_dense
    kpoints_sparse    = args.kpoints_sparse
    qpoints           = args.qpoints
    run_mode          = args.run_mode
    submit_job_system = args.submit_job_system
    # Running task types list
    run_mode_list = [
        "relax",                        # whether run relax.in or not
        "scffit",                       # whether run scf.fit.in or not
        "scf",                          # whether run scf.in or not
        "ph_no_split",                  # whether run run_ph_no_split.in or not
        "dyn0_flag",                    # whether only get *.dyn0 or not
        "ph_split_form_dyn0",           # whether run ph-split-from-dyn0.ph in the way1 or not
        "ph_split_set_startlast_q",     # whether run ph-split-startlastq.in in the way2 or not
        "q_non_irreducible_amount",     # non irreducible q amount you need input
        "merge",                        # whether run_merge or not
        "q2r",                          # whether run q2r.in or not
        "matdyn",                       # whether run matdyn.in or not
        "matdyn_dos",                   # whether run matdyn.dos.in or not
        "lambda"                        # whether run matdyn.dos.in or not
    ]

    if run_mode=="relax" and (work_path is not None):
        qe_sc = qe_workflow(
            input_file_path, 
            work_path, 
            kpoints_dense=kpoints_dense, 
            pressure=press, 
            run_mode=run_mode, 
            submit_job_system=submit_job_system)

    if run_mode=="scffit" and work_path is not None:
        # for root, dirs, files in os.walk(work_path):
        #     if 'relax.out' in files:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                kpoints_dense=kpoints_dense,
                run_mode=run_mode,
                submit_job_system=submit_job_system
                )
    
    if run_mode=="scf" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_workflow(relax_out_path, work_path, kpoints_sparse=kpoints_sparse)
                qe_sc.write_scf_in(root)
                if 'scf.in' in os.listdir(root):
                    qe_sc.slurmscf(root)
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system("sbatch slurmscf.sh")
                    logger.info("qe scf.fit is running")
                    os.chdir(cwd)
        
    if run_mode=="ph_no_split" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_workflow(relax_out_path, work_path, qpoints=qpoints)
                qe_sc.write_ph_no_split_in(root)
                if "ph_no_split.in" in os.listdir(root):
                    qe_sc.slurmph_no_split(root)
                    cwd = os.getcwd()
                    os.chdir(root)
                    id_info = os.popen("sbatch slurmph_no_split.sh").read()
                    if re.search(r"\d+", id_info) is not None:
                        id_num  = re.search(r"\d+", id_info).group()
                    os.chdir(cwd)
        while dyn0_flag:
            if os.path.exists(os.path.join(work_path, qe_sc.system_name+".dyn0")):
                os.system("sq")
                os.system("scancel {}".format(id_num))
                logger.info(f"The script detected the *.dyn0, so scancel the slurm job {id_num}")
                dyn0_flag = False
            else:
                time.sleep(5)
                dyn0_flag = True

    # run_ph_split_form_dyn0
    if run_mode=="ph_split_form_dyn0" and work_path is not None:
        files = os.listdir(work_path)
        dyn0_names = list(Path(work_path).glob("*.dyn0"))
        if 'relax.out' in files and 'scf.out' in files and dyn0_names:
            relax_out_path = os.path.join(work_path, 'relax.out')
            scf_out_path   = os.path.join(work_path, 'scf.out')
            if len(dyn0_names)==1:
                dyn0_path = dyn0_names[0]
            else:
                raise FileExistsError("Exist many *.dyn0 files or No *.dyn0")
            qe_sc = qe_workflow(relax_out_path, work_path, qpoints=qpoints)
            q_total_amount, q_non_irreducible_amount, q_coordinate_list = qe_sc.get_q_from_dyn0(dyn0_path)

        for i, q3 in enumerate(q_coordinate_list):
            split_ph_dir = os.path.join(work_path, str(i+1))
            if not os.path.exists(split_ph_dir):
                os.makedirs(split_ph_dir)
            qe_sc.write_split_ph_in_from_dyn0(split_ph_dir, q3)
            qe_sc.write_scf_fit_in(split_ph_dir)
            qe_sc.write_scf_in(split_ph_dir)
            
        for root, dirs, files in os.walk(work_path):
            if "split_ph.in" in files:
                qe_sc.slurmph_split_form_dyn0(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmph_split_form_dyn0.sh")
                os.chdir(cwd)

    # run_ph_split_set_startlast_q
    if run_mode=="ph_split_set_startlast_q" and work_path is not None:
        files = os.listdir(work_path)
        if 'relax.out' in files:
            relax_out_path = os.path.join(work_path, 'relax.out')
            # scf_out_path   = os.path.join(work_path, 'scf.out')
            qe_sc = qe_workflow(relax_out_path, work_path, qpoints=qpoints)
            # if  'scf.out' in files: 
            #     q_total_amount,    q_non_irreducible_amount, \
            #     q_coordinate_list, q_weight_list             = qe_sc.get_q_from_scfout(scf_out_path)
            if q_non_irreducible_amount != 0:
                for i in range(q_non_irreducible_amount):
                    qe_sc.write_split_ph_in_set_startlast_q(work_path, i+1, i+1)
        
        if work_path:
            split_ph_files = list(Path(work_path).glob("split_ph*.in"))
            if len(split_ph_files)==q_non_irreducible_amount:
                for split_ph_file in split_ph_files:
                    split_ph_name = str(split_ph_file).split(".")[0]
                    qe_sc.slurmph_split_set_startlast_q(work_path, split_ph_name)
                    cwd = os.getcwd()
                    os.chdir(work_path)
                    os.system("sbatch slurm_{}.sh".format(split_ph_name))
                    os.chdir(cwd)
                
    if run_mode=="merge" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if len(dirs) > 3 and "relax.out" in files:
                relax_out_path = os.path.join(root, 'relax.out')
                scf_out_path   = os.path.join(root, 'scf.out')
                qe_sc = qe_workflow(relax_out_path, work_path, qpoints=qpoints)
                dyn0_names = list(Path(work_path).glob("*.dyn0"))
                if len(dyn0_names)==1:
                    dyn0_path = dyn0_names[0]
                q_total_amount, q_non_irreducible_amount, q_coordinate_list = qe_sc.get_q_from_dyn0(dyn0_path)
                qe_sc.merge(root)

    if run_mode=="q2r" and work_path is not None:
        files = os.listdir(work_path)
        if 'relax.out' in files:
            relax_out_path = os.path.join(work_path, 'relax.out')
            qe_sc = qe_workflow(relax_out_path, work_path)
            qe_sc.write_q2r_in()
        if 'q2r.in' in os.listdir(work_path):
            qe_workflow.slurmq2r(work_path)
            cwd = os.getcwd()
            os.chdir(work_path)
            os.system("sbatch slurmq2r.sh")
            os.chdir(cwd)
    
    if run_mode=="matdyn" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_workflow(relax_out_path, work_path)
                qe_sc.write_matdyn_in()
            if 'matdyn.in' in os.listdir(root):
                qe_workflow.slurmmatgen(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmmatgen.sh")
                os.chdir(cwd)

    if run_mode=="matdyn_dos" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'matdyn.dos.in' in files:
                qe_workflow.slurmmatgen_dos(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmmatgen_dos.sh")
                os.chdir(cwd)

    if run_mode=="lambda" and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_workflow(relax_out_path, work_path, run_lambda=run_lambda)
            if 'lambda.in' in files:
                qe_workflow.slurmlambda(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmlambda.sh")
                os.chdir(cwd)
