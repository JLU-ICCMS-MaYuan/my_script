#!/usr/bin/env python3

from pathlib import Path

scripts_names = [
    'ares-phonon-toolkit', 
    'bader-toolkit',
    'calypso-toolkit', 
    'lobster-toolkit', 

    'deepkit-dpgen-toolkit',
    'ELK-toolkit',
    'EPW-toolkit',

    'mytoolkit/convexhull', 
    'mytoolkit/energy_curve_with_pressure', 
    'mytoolkit/files_tidy',
    'mytoolkit/job_tasks_submit',
    'mytoolkit/pieces',
    'mytoolkit/qe_toolkit',
    'mytoolkit/structures_characters',
    'mytoolkit/superconductivity',
    'mytoolkit/vasp_toolkit',
    'mytoolkit/vasp_toolkit/vasp-MD-scripts',
    'mytoolkit/wien2K',

    'uspex-toolkit', 
    ]


with open(Path('~/.myenvs').expanduser(), 'w') as f:
    for dirname in scripts_names:
        abs_path = Path(dirname).absolute()
        f.write(f'export PATH="{abs_path}:$PATH"\n')
        
