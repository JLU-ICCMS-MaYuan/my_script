from pathlib import Path

from argparse import ArgumentParser

from mytoolkit.convert import convert
from mytoolkit.config import config
from mytoolkit.kmesh import create_kmesh

def format_convert(args: ArgumentParser):
    config_d = config(args).read_config()
    convert(
        input_file_path=config_d.get('input_file_path', None),
        work_path=config_d.get('work_path', './'),
        dst_format=config_d.get('dst_format', None),
    )

def kmesh(args: ArgumentParser):
    config_d = config(args).read_config()
    input(config_d)
    kpoints_path = Path(config_d.get('work_path', './')).joinpath("KPOINTS") 
    create_kmesh(
        kresolution=config_d.get('kspacing', 0.4),
        input_file_path=config_d.get('input_file_path', None),
        output_kpoints=kpoints_path,
    )