#!/usr/bin/env python3
"""
从指定的 IT 目录的 cam 文件中提取结构，并复制对应的 .res 文件。
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="提取 cam 指定行范围内的结构，可选只保留三元结构。"
    )
    parser.add_argument(
        "it",
        help="IT 编号，既可以写成 IT5 也可以写成 5。",
    )
    parser.add_argument(
        "-b",
        "--begin",
        type=int,
        required=True,
        help="开始行号（从 1 起算，包含该行）。",
    )
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        required=True,
        help="结束行号（从 1 起算，包含该行）。",
    )
    parser.add_argument(
        "-t",
        "--ternary-only",
        action="store_true",
        help="只保留三元（含三个元素）结构。",
    )
    return parser.parse_args()


def normalize_it(raw_it: str) -> str:
    if not raw_it:
        raise ValueError("IT 编号不能为空")

    it = raw_it.strip()
    prefix = "IT"
    if it.lower().startswith(prefix.lower()):
        suffix = it[len(prefix) :]
    else:
        suffix = it

    if not suffix.isdigit():
        raise ValueError(f"无法解析 IT 编号：{raw_it}")

    return f"{prefix}{int(suffix)}"


def count_unique_elements(formula: str) -> int:
    pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
    return len({match.group(1) for match in pattern.finditer(formula)})


def main() -> int:
    args = parse_args()

    if args.begin < 1 or args.end < 1:
        print("开始和结束行号必须大于等于 1。", file=sys.stderr)
        return 1
    if args.begin > args.end:
        print("开始行号不能大于结束行号。", file=sys.stderr)
        return 1

    try:
        it_name = normalize_it(args.it)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    base_dir = Path(__file__).resolve().parent
    it_dir = base_dir / it_name
    if not it_dir.is_dir():
        print(f"未找到目录：{it_dir}", file=sys.stderr)
        return 1

    cam_path = it_dir / "cam"
    if not cam_path.is_file():
        print(f"未找到 cam 文件：{cam_path}", file=sys.stderr)
        return 1

    lines = cam_path.read_text().splitlines()
    total_lines = len(lines)
    if args.end > total_lines:
        print(
            f"结束行号 {args.end} 超出 cam 总行数 {total_lines}。",
            file=sys.stderr,
        )
        return 1

    save_dir = base_dir / f"save_structs_{it_name}"
    save_dir.mkdir(parents=True, exist_ok=True)

    copied = []
    missing = []
    written_lines = []

    for idx, raw_line in enumerate(lines, start=1):
        if idx < args.begin or idx > args.end:
            continue

        line = raw_line.rstrip()
        if not line:
            continue

        parts = line.split()
        if not parts:
            continue

        struct_id = parts[0]
        composition = parts[8] if len(parts) > 8 else ""

        if args.ternary_only:
            if not composition:
                continue
            if count_unique_elements(composition) != 3:
                continue

        res_path = it_dir / f"{struct_id}.res"
        if not res_path.is_file():
            missing.append((idx, struct_id))
            continue

        dest_path = save_dir / res_path.name
        shutil.copy2(res_path, dest_path)
        copied.append((idx, struct_id))
        written_lines.append(f"{idx}\t{line}")

    if written_lines:
        summary_path = save_dir / "selected_lines.tsv"
        summary_path.write_text("\n".join(written_lines) + "\n")

    print(f"输出目录：{save_dir}")
    print(f"成功复制 {len(copied)} 个结构。")
    if missing:
        print("以下结构未找到对应的 .res 文件：", file=sys.stderr)
        for idx, name in missing:
            print(f"  行 {idx}: {name}", file=sys.stderr)

    if args.ternary_only and not copied:
        print("在指定范围内未找到符合三元条件的结构。", file=sys.stderr)

    return 0 if copied or not args.ternary_only else 1


if __name__ == "__main__":
    sys.exit(main())

