from argparse import ArgumentParser, RawTextHelpFormatter
from asyncio import subprocess

from mytoolkit.tool_run import format_convert, kmesh

def set_more_args(parser: ArgumentParser):
    # 指明输入文件所在的位置
    parser.add_argument(
        '-i',
        '--input-file-path',
        type=str,
        default=None,
        dest='input_file_path',
        help="please tell me your inputfile path\n"
            '请一定注意, 当使用kmesh子功能时, 输入文件必须是POSCAR\n'
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
    parser_formatconvert = subparsers.add_parser("convert", formatter_class=RawTextHelpFormatter)
    parser_formatconvert.add_argument(
        '-m',
        '--more-args-about-what-destational-format-you-want',
        type=str,
        nargs="+",
        dest='more_args',
        help='你的目标格式\n'
            '例如: convert -m dst_format=cif\n'
            '例如: convert -m dst_format=vasp\n'
            '例如: convert -m dst_format=wien2k\n'
    )
    parser_formatconvert.set_defaults(tool_flow=format_convert)

    # 生成KPOINTS
    parser_kmesh = subparsers.add_parser("kmesh", formatter_class=RawTextHelpFormatter)
    parser_kmesh.add_argument(
        '-m',
        '--more-args-about-kspacing-meshing-density',
        type=str,
        nargs="+",
        dest='more_args',
        help='请指定你的K点网格密度\n'
            '例如: kmesh -m kspacing=0.1\n'
            '请一定注意输入文件必须是POSCAR\n'
    )
    parser_kmesh.set_defaults(tool_flow=kmesh)

    # 绘制 压制相图
    parser_plot_phase_transition = subparsers.add_parser("plot_phase_transition", formatter_class=RawTextHelpFormatter)
    parser_plot_phase_transition.add_argument(
        '-m',
        '--more-args-about-plot-phase-transition',
        type=str,
        nargs="+",
        dest='more_args',
        help='请指定\n'
            '例如: kmesh -m kspacing=0.1\n'
            '请一定注意输入文件必须是POSCAR\n'
    )

    args = parser.parse_args()
    return args