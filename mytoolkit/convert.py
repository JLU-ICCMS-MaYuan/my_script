from pathlib import Path

from ase.io import read, write


def convert(
    input_file_path,
    work_path,
    dst_format,
    ):

    if input_file_path is None or dst_format is None:
        raise ValueError("你必须指定输入文件路径和输出文件格式")

    input_file_path = Path(input_file_path)
    work_path       = Path(work_path)
    struct          = read(input_file_path)
    formula         = struct.get_chemical_formula("metal")

    if dst_format == "cif":
        dst_file = work_path.joinpath(formula+".cif")
        struct.write(dst_file, format="cif")
    elif dst_format == "vasp":
        dst_file = work_path.joinpath(formula+".vasp")
        struct.write(dst_file, format="vasp")