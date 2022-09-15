from argparse import ArgumentParser
from convert import convert
from config import config


def format_convert(args: ArgumentParser):
    config_d = config(args).read_config()
    convert(
        input_file_path=config_d.get('input_file_path', None),
        work_path=config_d.get('work_path', './'),
        dst_format=config_d.get('dst_format', None),
    )