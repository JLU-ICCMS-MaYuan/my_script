#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
#!/work/home/mayuan/miniconda3/envs/cage/bin/python3

"""
# check tasks number
j=0;x=1; for i in `squeue | awk '{print $1}'`; do  let j+=x; done; echo $j
pbs
    qe_main.py -i ./test/POSCAR -w ./test/out -kd 16 16 16 -p 150 -mode relax -j pbs 
    qe_main.py -w ./test/out/150.0/  -kd 16 16 16 -mode scffit -j pbs
    qe_main.py -w ./test/out/150.0/  -kd 8 8 8 -mode scf -j pbs

    # no split method
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j pbs -mode ph_no_split 
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j pbs -mode ph_no_split -dyn0

    # spilit method 1
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j pbs -mode ph_no_split -dyn0
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j pbs -mode ph_split_form_dyn0 
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -mode merge

    # split method 2
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j pbs -mode ph_split_set_startlast_q -q-no-irr-num 8

    qe_main.py -w ./test/out/150.0/ -j pbs -mode q2r 
    qe_main.py -w ./test/out/150.0/ -j pbs -mode matdyn
    qe_main.py -w ./test/out/150.0/ -j pbs -mode matdyn_dos
    qe_main.py -w ./test/out/150.0/ -j pbs -mode lambda -q 4 4 4

slurm
    qe_main.py -i ./test/POSCAR -w ./test/out -kd 16 16 16 -p 150 -mode relax -j slurm
    qe_main.py -w ./test/out/150.0/  -kd 16 16 16 -mode scffit -j slurm
    qe_main.py -w ./test/out/150.0/  -ks 8 8 8 -mode scf -j slurm

    # no split method
    qe_main.py -w ./test/out/150.0/ -mode ph_no_split -q 8 8 8 -j slurm
    qe_main.py -w ./test/out/150.0/ -mode ph_no_split -q 8 8 8 -j slurm -dyn0

    # spilit method 1
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j slurm -mode ph_no_split  -dyn0
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j slurm -mode ph_split_form_dyn0 
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -mode merge

    # split method 2
    qe_main.py -w ./test/out/150.0/ -q 4 4 4 -j slurm -mode ph_split_set_startlast_q -q-no-irr-num 8

    qe_main.py -w ./test/out/150.0/ -j slurm -mode q2r 
    qe_main.py -w ./test/out/150.0/ -j slurm -mode matdyn
    qe_main.py -w ./test/out/150.0/ -j slurm -mode matdyn_dos
    qe_main.py -w ./test/out/150.0/ -j slurm -mode lambda -q 4 4 4
"""
import os
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
        '-q-no-irr-num',
        '--q-non-irreducible-amount',
        action='store',
        type=int,
        default=0,
        dest='q_non_irreducible_amount',
        help="non irreducible q amount you need input",
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
    parser.add_argument(
        '-dyn0',
        '--for-get-dyn0',
        action='store_true',
        default=False,
        dest='dyn0_flag',
        help="whether only get *.dyn0 or not",
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

    dyn0_flag                = args.dyn0_flag
    q_non_irreducible_amount = args.q_non_irreducible_amount

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
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                kpoints_sparse=kpoints_sparse,
                run_mode=run_mode,
                submit_job_system=submit_job_system
                )
        
    if run_mode=="ph_no_split" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
                submit_job_system=submit_job_system,
                dyn0_flag=dyn0_flag
                )
    # run_ph_split_form_dyn0
    if run_mode=="ph_split_form_dyn0" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
                submit_job_system=submit_job_system,
                )

    # run_ph_split_set_startlast_q
    if run_mode=="ph_split_set_startlast_q" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
                submit_job_system=submit_job_system,
                q_non_irreducible_amount=q_non_irreducible_amount
                )

    if run_mode=="merge" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
            )
    
    if run_mode=="q2r" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                run_mode=run_mode,
                submit_job_system=submit_job_system,
            )
    
    if run_mode=="matdyn" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
                submit_job_system=submit_job_system,
            )

    if run_mode=="matdyn_dos" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                run_mode=run_mode,
                submit_job_system=submit_job_system,
            )

    if run_mode=="lambda" and work_path is not None:
        relax_out = list(Path(work_path).glob("relax.out"))
        if relax_out:
            relax_out_path = str(relax_out[0].absolute())
            qe_sc = qe_workflow(
                relax_out_path, 
                work_path, 
                qpoints=qpoints,
                run_mode=run_mode,
                submit_job_system=submit_job_system,
            )
