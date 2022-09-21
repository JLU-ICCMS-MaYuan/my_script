from argparse import ArgumentParser
from asyncio import subprocess

from tool_run import format_convert

def set_more_args(parser: ArgumentParser):
    # 指明输入文件所在的位置
    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your inputfile path",
    )
    # 指明工作目录
    parser.add_argument(
        '-w',
        '-work-path',
        type=str,
        default=None,
        dest='work_path',
        help="please tell me your calculated directory, I will put all file there",
    )


    subparsers = parser.add_subparsers()


    # 格式转化
    parser_formatconvert = subparsers.add_parser("convert")
    parser_formatconvert.add_argument(
        '-m',
        '--more-args-about-what-destational-format-you-want',
        type=str,
        dest='more_args',
        help='你的目标格式',
    )
    parser_formatconvert.set_defaults(tool_flow=format_convert)



    args = parser.parse_args()
    return args