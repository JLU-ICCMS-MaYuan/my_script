#!/usr/bin/env python3

from pathlib import Path

scripts_names = [
    'ares-phonon-toolkit', 
    'bader-toolkit',
    'calypso-toolkit', 
    'lobster-toolkit', 

    'deepkit-dpgen-toolkit',
    'elk-toolkit',
    'epw/epw-toolkit',
    'lammps-toolkit',

    'mytoolkit/convexhull', 
    'mytoolkit/energy_curve_with_pressure', 
    'mytoolkit/files_tidy',
    'mytoolkit/job_tasks_submit/slurm',
    'mytoolkit/job_tasks_submit/pbs',
    'mytoolkit/job_tasks_submit/lsf',
    'mytoolkit/pieces',
    'mytoolkit/qe_toolkit',
    'mytoolkit/structures_characters',
    'mytoolkit/superconductivity',
    'mytoolkit/vasp_toolkit',
    'mytoolkit/vasp_toolkit/vasp-MD-scripts',
    'mytoolkit/wien2K',
    'mytoolkit/many_jobsubmitting_system/slurm',
    'mytoolkit/many_jobsubmitting_system/pbs',
    'mytoolkit/many_jobsubmitting_system/lsf',

    'uspex-toolkit', 
    'magus-toolkit',
    'sscha-toolkit',
    'phonopy-toolkit',
    ]


with open(Path('~/.myenvs').expanduser(), 'w') as f:
    for dirname in scripts_names:
        abs_path = Path(dirname).absolute()
        f.write(f'export PATH="{abs_path}:$PATH"\n')
        
