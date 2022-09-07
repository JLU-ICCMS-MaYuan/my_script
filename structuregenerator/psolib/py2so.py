import argparse
import shutil
import sys
from pathlib import Path

import setuptools
from Cython.Build import cythonize

# python py2so.py [src] [build]


def copyfile(src, root, dst):
    """Copy src to dst, keep relative path to root"""
    f = Path(src).relative_to(root)
    new_f = Path(dst).joinpath(f)
    new_f.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, new_f)
    print(f"Copying {f} -> {new_f}")


def get_py(
    root_dir: Path,
    src: str,
    build: str,
    include_patterns=[],
    exclude_patterns=[],
):
    src_dir = Path(root_dir).joinpath(src)
    build_dir = Path(root_dir).joinpath(build)
    for f in src_dir.rglob('*'):
        if any([f.match(pat) for pat in include_patterns]):
            copyfile(f, root_dir, build_dir)
        elif any([f.match(pat) for pat in exclude_patterns]):
            pass
        elif f.suffix in ['.py', '.pyx']:
            yield str(f)
        else:
            print(f"Ignore {f}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser.add_argument('src', nargs='+', help="Source files or directories")
    parser.add_argument('build', help="Build directory")
    parser.add_argument(
        '--include-patterns',
        default=['main.py', '__init__.py', '*.so'],
        nargs='+',
        help="Do Not compile and directly copy to build directory",
    )
    parser.add_argument(
        '--exclude-patterns',
        default=['__pycache__/*', '*.c', '*.o'],
        nargs='+',
        help="Ignored files patterns",
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    root_dir = Path('.')
    include_patterns = ['main.py', '__init__.py', '*.so'] + args.include_patterns
    exclude_patterns = list(
        set(['__pycache__/*', '*.c', '*.o'] + args.exclude_patterns)
    )
    py_list = []
    for src in args.src:
        if Path(src).is_dir():
            py_list += list(
                get_py(
                    root_dir,
                    src,
                    args.build,
                    include_patterns,
                    exclude_patterns,
                )
            )
        elif any([Path(src).match(pat) for pat in include_patterns]):
            copyfile(src, root_dir, args.build)
        elif Path(src).is_file() and Path(src).suffix == '.py':
            py_list.append(src)
        else:
            print(f"Ignore {src}")

    build_tmp = Path(args.build) / 'tmp'

    setuptools.setup(
        ext_modules=cythonize(py_list),
        script_args=['build_ext', '-b', f'{args.build}', '-t', f'{build_tmp}'],
    )

    if Path(build_tmp).exists():
        shutil.rmtree(build_tmp)
