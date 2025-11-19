#!/usr/bin/env python3

"""
批量将多个 VASP 目录（包含 OUTCAR）转换为 deepmd-format 原始/NPY 数据。
示例：
    python convert_vasp_dirs.py ../run1 ../run2 --raw-dir 1.traindata --npy-dir 1.traindata --set-size 5
"""

import argparse
import glob
import os
from pathlib import Path

from dpdata import LabeledSystem, MultiSystems


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="批量把 VASP OUTCAR 转换为 deepmd raw 与 npy 数据"
    )
    parser.add_argument("inputs", nargs="*", help="包含 OUTCAR 的目录，可指定多个")
    parser.add_argument(
        "--set-size",
        type=int,
        default=5,
        help="to_deepmd_npy 的 set_size 参数（默认：5）",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="1.traindata",
        help="输出 raw/npy(set.*) 的目录（默认：1.traindata）",
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="保存 OUTCAR 路径列表的文件（每行一个路径，允许空行与通配符）",
    )
    return parser.parse_args()


def _collect_from_dirs(directories: list[str]) -> list[Path]:
    """从目录列表中收集 OUTCAR 路径。"""
    results: list[Path] = []
    for directory in directories:
        outcar_path = Path(directory).expanduser().resolve() / "OUTCAR"
        results.append(outcar_path)
    return results


def _collect_from_file(file_path: str) -> list[Path]:
    """从包含 OUTCAR 路径的文件中读取路径，支持通配符。"""
    path_obj = Path(file_path).expanduser().resolve()
    if not path_obj.exists():
        raise SystemExit(f"[ERROR] 指定的文件 {path_obj} 不存在。")
    lines = path_obj.read_text().splitlines()
    collected: list[Path] = []
    for idx, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        expanded_pattern = os.path.expandvars(os.path.expanduser(line))
        matches = glob.glob(expanded_pattern, recursive=True)
        if not matches:
            print(f"[WARN] 行 {idx}: {line} 未匹配到任何路径，跳过。")
            continue
        for match in matches:
            match_path = Path(match).expanduser().resolve()
            if match_path.is_dir():
                raise SystemExit(
                    f"[ERROR] 行 {idx}: {match_path} 是目录，预期为 OUTCAR 文件。"
                )
            collected.append(match_path)
    return collected


def main() -> None:
    args = parse_args()

    if not args.inputs and not args.file:
        raise SystemExit("必须提供至少一个目录或指定 --file。")

    outcar_candidates: list[Path] = []
    if args.inputs:
        outcar_candidates.extend(_collect_from_dirs(args.inputs))
    if args.file:
        outcar_candidates.extend(_collect_from_file(args.file))

    # 去重并保持顺序
    unique_outcars: list[Path] = []
    seen: set[str] = set()
    for candidate in outcar_candidates:
        key = str(candidate)
        if key in seen:
            continue
        unique_outcars.append(candidate)
        seen.add(key)

    total_datas = MultiSystems()
    for outcar_path in unique_outcars:
        print(f"---------- convert {outcar_path} -> dpdata -----------")
        if not outcar_path.exists():
            print(f"[WARN] 找不到 {outcar_path}，跳过。")
            continue
        if outcar_path.is_dir():
            raise SystemExit(f"[ERROR] {outcar_path} 是目录，预期为 OUTCAR 文件。")
        try:
            ls = LabeledSystem(str(outcar_path), fmt="outcar")
            total_datas.append(ls)
        except ValueError as exc:
            print(f"[ERROR] 解析 {outcar_path} 失败：{exc}")

    if len(total_datas.systems) == 0:
        raise SystemExit("未成功解析任何 OUTCAR，退出。")

    out_dir = Path(args.output).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    total_datas.to_deepmd_raw(str(out_dir))
    total_datas.to_deepmd_npy(str(out_dir), set_size=args.set_size)
    print(f"完成：raw/npy 输出到 {out_dir}，set_size={args.set_size}")


if __name__ == "__main__":
    main()
