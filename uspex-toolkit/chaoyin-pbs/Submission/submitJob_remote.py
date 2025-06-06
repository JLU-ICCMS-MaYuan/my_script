

import argparse
import os
import re

from subprocess import check_output
from io import open

_author_ = 'etikhonov'


def submitJob_remote(workingDir, index, commandExecutable):
    """
    This routine is to submit job to remote cluster

    One needs to do a little edit based on your own case.

    Step 1: to prepare the job script which is required by your supercomputer
    Step 2: to submit the job with the command like qsub, bsub, llsubmit, .etc.
    Step 3: to get the jobID from the screen message

    :param workingDir: working directory on remote machine
    :param index: index of the structure.
    :param commandExecutable: command executable for current step of optimization
    :return:
    """

    # Step 1
    # Specify the PATH to put your calculation folder
    Home = '/home/etikhonov' # 'pwd' of your home directory of your remote machine
    Address = 'rurik'  # your target server: ssh alias or username@address
    Path = Home + '/' + workingDir + '/CalcFold' + -str(index) # Just keep it
    run_content = ''
    run_content += '#!/bin/sh\n'
    run_content += '#PBS -N    opt' + str(index) + '\n'
    run_content += '#PBS -q    liuhy\n'
    run_content += '#PBS -l    nodes=1:ppn=28\n'
    run_content += '#PBS -j    oe\n'
    run_content += '#PBS -V\n'
    run_content += 'source /public/home/mayuan/intel/oneapi/setvars.sh --force\n'
    run_content += 'cd $PBS_O_WORKDIR'

    run_content += 'ulimit -s unlimited\n'
    run_content += 'export I_MPI_ADJUST_REDUCE=3\n'
    run_content += 'export MPIR_CVAR_COLL_ALIAS_CHECK=0\n'
    run_content += 'export I_MPI_FABRICS=shm\n'
    run_content += 'export MKL_DEBUG_CPU_TYPE=5\n'

    run_content += 'cd '+ Path + '\n'
    run_content += commandExecutable + '\n'

    with open('myrun', 'w') as fp:
        fp.write(run_content)

    # Create the remote directory
    # Please change the ssh/scp command if necessary.
    try:
        os.system('ssh -i ~/.ssh/id_rsa ' + Address + ' mkdir -p ' + Path)
    except:
        pass

    # Copy calculation files
    # add private key -i ~/.ssh/id_rsa if necessary
    os.system('scp POSCAR   ' + Address + ':' + Path)
    os.system('scp INCAR    ' + Address + ':' + Path)
    os.system('scp POTCAR   ' + Address + ':' + Path)
    os.system('scp KPOINTS  ' + Address + ':' + Path)
    os.system('scp myrun ' + Address + ':' + Path)

    # Step 2
    # Run command
    output = str(check_output('ssh -i ~/.ssh/id_rsa ' + Address + ' qsub ' + Path + '/myrun', shell=True))

    # Step 3
    # Here we parse job ID from the output of previous command
    jobNumber = int(re.findall(r'\d+', output)[0])
    return jobNumber


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='index', type=int)
    parser.add_argument('-c', dest='commandExecutable', type=str)
    parser.add_argument('-f', dest='workingDir', type=str)
    args = parser.parse_args()

    jobNumber = submitJob_remote(workingDir=args.workingDir, index=args.index, commandExecutable=args.commandExecutable)
    print('<CALLRESULT>')
    print((int(jobNumber)))
