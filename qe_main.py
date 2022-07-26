#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

"""
# check tasks number
j=0;x=1; for i in `squeue | awk '{print $1}'`; do  let j+=x; done; echo $j

relax:24 24 24
    qe_main.py -i ./POSCAR -w ./out -relax -kd 24 24 24 -p 200 

    qe_main.py -w ./ -scf -ks 6 6 6

    qe_main.py -w ./ -scffit -kd 12 12 12 

    qe_main.py -w ./ -ph-split -kd 24 24 24 

    qe_main.py -w ./ -m



    qe_main.py -w ./ -ph-nosplit -q 4 4 4

    qe_main.py -w ./ -ph-nosplit -q 8 8 8 -dyn0

"""

import os
import re
import shutil
import time
import logging
from argparse import ArgumentParser

from qeSuperconductTc import qe_superconduct_workflow as qe_sc_workflow

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
        '-relax',
        '--run-relax',
        action='store_true',
        default=False,
        dest='run_relax',
        help="whether run relax.in or not",
    )
    parser.add_argument(
        '-scffit',
        '--run-scf-fit',
        action='store_true',
        default=False,
        dest='run_scfFit',
        help="whether run scf.fit.in or not",
    )
    parser.add_argument(
        '-scf',
        '--run-scf',
        action='store_true',
        default=False,
        dest='run_scf',
        help="whether run scf.in or not",
    )
    parser.add_argument(
        '-ph-nosplit',
        '--run-ph-no-split',
        action='store_true',
        default=False,
        dest='run_ph_no_split',
        help="whether run run_ph_no_split.in or not",
    )
    parser.add_argument(
        '-dyn0',
        '--for-get-dyn0',
        action='store_true',
        default=False,
        dest='dyn0_flag',
        help="whether only get *.dyn0 or not",
    )
    parser.add_argument(
        '-ph-split',
        '--run-ph-split',
        action='store_true',
        default=False,
        dest='run_ph_split',
        help="whether run ph.in or not",
    )
    parser.add_argument(
        '-m',
        '--run-merge',
        action='store_true',
        default=False,
        dest='run_merge',
        help="whether run_merge or not",
    )
    parser.add_argument(
        '-q2r',
        '--run-q2r',
        action='store_true',
        default=False,
        dest='run_q2r',
        help="whether run q2r.in or not",
    )
    parser.add_argument(
        '-matdyn',
        '--run-matdyn',
        action='store_true',
        default=False,
        dest='run_matdyn',
        help="whether run matdyn.in or not",
    )
    parser.add_argument(
        '-dos',
        '--run-matdyn-dos',
        action='store_true',
        default=False,
        dest='run_matdyn_dos',
        help="whether run matdyn.dos.in or not",
    )
    parser.add_argument(
        '-lambda',
        '--run-lambda',
        action='store_true',
        default=False,
        dest='run_lambda',
        help="whether run matdyn.dos.in or not",
    )
    args = parser.parse_args()
    input_file_path = args.input_file_path
    press           = args.press
    work_path       = args.work_path

    kpoints_dense   = args.kpoints_dense
    kpoints_sparse  = args.kpoints_sparse
    qpoints         = args.qpoints

    run_relax       = args.run_relax
    run_scfFit      = args.run_scfFit
    run_scf         = args.run_scf
    
    run_ph_no_split = args.run_ph_no_split
    dyn0_flag       = args.dyn0_flag

    run_ph_split    = args.run_ph_split
    run_merge       = args.run_merge
    run_q2r         = args.run_q2r
    run_matdyn      = args.run_matdyn
    run_matdyn_dos  = args.run_matdyn_dos      
    run_lambda      = args.run_lambda

    if run_relax and work_path is not None:
        qe_sc = qe_sc_workflow(input_file_path, work_path, kpoints_dense=kpoints_dense, pressure=press, run_relax=run_relax)
        for root, dirs, files in os.walk(work_path):
            if 'relax.in' in files:
                qe_sc_workflow.slurmrelax(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmrelax.sh")
                logger.info("qe relax is running")
                os.chdir(cwd)

    if run_scfFit and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, kpoints_dense=kpoints_dense, run_scfFit=run_scfFit)
                if 'scf.fit.in' in os.listdir(root):
                    qe_sc.slurmscfFit(root)
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system("sbatch slurmscfFit.sh")
                    logger.info("qe scf.fit is running")
                    os.chdir(cwd)

    if run_scf and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                if kpoints_sparse == [1, 1, 1]:
                    qe_sc = qe_sc_workflow(relax_out_path, work_path, kpoints_dense=kpoints_dense, run_scf=run_scf)
                else:
                    qe_sc = qe_sc_workflow(relax_out_path, work_path, kpoints_sparse=kpoints_sparse, run_scf=run_scf)
                if 'scf.in' in os.listdir(root):
                    qe_sc.slurmscf(root)
                    cwd = os.getcwd()
                    os.chdir(root)
                    os.system("sbatch slurmscf.sh")
                    logger.info("qe scf.fit is running")
                    os.chdir(cwd)
           
    if run_ph_no_split and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, qpoints=qpoints, run_ph_no_split=run_ph_no_split)
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

    if run_ph_split and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files and 'scf.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                scf_out_path   = os.path.join(root, 'scf.out')
                qe_sc = qe_sc_workflow(
                                    relax_out_path, work_path, 
                                    kpoints_dense=kpoints_dense, 
                                    run_ph_split=run_ph_split, 
                                    scf_out_path=scf_out_path
                                    )
        for root, dirs, files in os.walk(work_path):
            if "split_ph.in" in files:
                qe_sc.slurmph_split(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmph_split.sh")
                os.chdir(cwd)
        
    if run_ph_split and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(
                                    relax_out_path, work_path, 
                                    kpoints_dense=kpoints_dense, 
                                    run_ph_split=run_ph_split,
                                    )

    if run_merge and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if len(dirs) > 3 and "relax.out" in files and 'scf.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                scf_out_path   = os.path.join(root, 'scf.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, qpoints=qpoints)
                qe_sc.get_q(scf_out_path)
                qe_sc.merge(root)

    if run_q2r and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, run_q2r=run_q2r)
            if 'q2r.in' in os.listdir(root):
                qe_sc_workflow.slurmq2r(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmq2r.sh")
                os.chdir(cwd)
    
    if run_matdyn and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, run_matdyn=run_matdyn)
            if 'matdyn.in' in os.listdir(root):
                qe_sc_workflow.slurmmatgen(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmmatgen.sh")
                os.chdir(cwd)

    if run_matdyn_dos and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'matdyn.dos.in' in files:
                qe_sc_workflow.slurmmatgen_dos(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmmatgen_dos.sh")
                os.chdir(cwd)

    if run_lambda and work_path is not None:
        for root, dirs, files in os.walk(work_path):
            if 'relax.out' in files:
                relax_out_path = os.path.join(root, 'relax.out')
                qe_sc = qe_sc_workflow(relax_out_path, work_path, run_lambda=run_lambda)
            if 'lambda.in' in files:
                qe_sc_workflow.slurmlambda(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmlambda.sh")
                os.chdir(cwd)