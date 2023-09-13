from __future__ import with_statement
from __future__ import absolute_import
from subprocess import check_output
import re
import sys
from io import open



def submitJob_local(index, commandExecutable):
    """
    This routine is to submit job locally
    One needs to do a little edit based on your own case.

    Step 1: to prepare the job script which is required by your supercomputer
    Step 2: to submit the job with the command like qsub, bsub, llsubmit, .etc.
    Step 3: to get the jobID from the screen message
    :return: job ID
    """

    RUN_FILENAME = 'myrun'
    JOB_NAME = 'optrun-{}'.format(index)

    # Step 1
    myrun_content = '''#!/bin/sh
#!/bin/bash
#SBATCH  --job-name={}
#SBATCH  --output=opt.out
#SBATCH  --error=opt.err                                  
#SBATCH  --partition=lhy
#SBATCH  --nodes=1                                                       
#SBATCH  --ntasks=48                                                
#SBATCH  --ntasks-per-node=48                                            
#SBATCH  --cpus-per-task=1   
 
source /work/env/intel2018
  
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0


{}
'''.format(JOB_NAME, commandExecutable)


    with open(RUN_FILENAME, 'w') as fp:
        fp.write(myrun_content)

    # Step 2
    # It will output some message on the screen like '2350873.nano.cfn.bnl.local'
    output = check_output('sbatch {}'.format(RUN_FILENAME), shell=True, universal_newlines=True)


    # Step 3
    # Here we parse job ID from the output of previous command
    jobNumber = int(output.split(' ')[3])
    print(str(jobNumber))
    return jobNumber


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='index', type=int)
    parser.add_argument('-c', dest='commandExecutable', type=str)
    args = parser.parse_args()

    jobNumber = submitJob_local(index=args.index, commandExecutable=args.commandExecutable)
    print('<CALLRESULT>')
    print(int(jobNumber))
