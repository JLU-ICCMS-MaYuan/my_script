#!/usr/bin/env python3

"""
批量将多个 VASP 目录（包含 OUTCAR）转换为 deepmd-format 原始/NPY 数据。
示例：
    python convert_vasp_dirs.py ../run1 ../run2 --raw-dir 1.traindata --npy-dir 1.traindata --set-size 5
"""

import argparse
from pathlib import Path

from dpdata import LabeledSystem, MultiSystems


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="批量把 VASP OUTCAR 转换为 deepmd raw 与 npy 数据"
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="包含 OUTCAR 的目录，可指定多个",
    )
    parser.add_argument(
        "--set-size",
        type=int,
        default=5,
        help="to_deepmd_npy 的 set_size 参数（默认：5）",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    total_datas = MultiSystems()
    for directory in args.inputs:
        outcar_path = Path(directory).expanduser().resolve() / "OUTCAR"
        print(f"---------- convert {outcar_path} -> dpdata -----------")
        if not outcar_path.exists():
            print(f"[WARN] 找不到 {outcar_path}，跳过。")
            continue
        try:
            ls = LabeledSystem(str(outcar_path), fmt="outcar")
            total_datas.append(ls)
        except ValueError as exc:
            print(f"[ERROR] 解析 {outcar_path} 失败：{exc}")

    if len(total_datas.systems) == 0:
        raise SystemExit("未成功解析任何 OUTCAR，退出。")

    out_dir = Path("1.traindata")
    total_datas.to_deepmd_raw(str(out_dir))
    total_datas.to_deepmd_npy(str(out_dir), set_size=args.set_size)
    print(f"完成：raw/npy 输出到 {out_dir}，set_size={args.set_size}")


if __name__ == "__main__":
    main()

