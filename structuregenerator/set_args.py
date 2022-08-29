from argparse import ArgumentParser
from ast import parse

from generator_methods import generator_methods

def set_more_args(parser: ArgumentParser):

    parser.add_argument(
        "-w",
        "--work-path",
        action="store",
        dest="work_path",
        type=str,
        help="specify the output path for news strutures"
    )

    parser.add_argument(
        "-i",
        "--input-file-path",
        action="store",
        default=str,
        dest="input_file_path", 
        help="specify the input file path",
    )

    subparsers = parser.add_subparsers(help="subparsers")
    parser_genmethod = subparsers.add_parser("method")
    parser_genmethod.set_defaults(generator=generator_methods)
    parser_genmethod.add_argument(
        '-m',
        '--more-argments-about-generate-struct-method',
        type=str,
        dest='more_args',
        nargs='+',
        help="输入更多关于结构弛豫的参数"
    )

    args = parser.parse_args()

    return args