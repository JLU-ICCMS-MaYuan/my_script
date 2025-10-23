#!/usr/bin/env python3
"""统计指定迭代目录中各个元素配比出现的次数。"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Iterable, List


def parse_iteration_label(label: str) -> str:
    value = label.strip()
    if not value:
        raise ValueError("迭代编号不能为空")
    upper = value.upper()
    if upper.startswith("IT"):
        return upper
    return f"IT{upper}"


def gcd_many(values: Iterable[int]) -> int:
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration:
        return 1
    result = abs(result)
    for value in iterator:
        result = math.gcd(result, abs(value))
        if result == 1:
            return 1
    return max(result, 1)


def parse_composition(path: Path) -> OrderedDict[str, int]:
    counts: Counter[str] = Counter()
    sfac_order: List[str] = []
    in_atoms = False
    with path.open("r", encoding="ascii", errors="ignore") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            upper = line.upper()
            if upper.startswith("TITL") or upper.startswith("CELL") or upper.startswith("LATT"):
                continue
            if upper.startswith("SFAC"):
                parts = raw_line.split()
                for token in parts[1:]:
                    token_clean = token.strip()
                    if not token_clean:
                        continue
                    sfac_order.append(token_clean)
                in_atoms = True
                continue
            if upper.startswith("UNIT") or upper.startswith("REM"):
                continue
            if upper == "END":
                break
            if not in_atoms:
                continue
            parts = raw_line.split()
            if not parts:
                continue
            element = parts[0].strip()
            if not element or element.upper() == "END":
                continue
            counts[element] += 1
    ordered_counts: OrderedDict[str, int] = OrderedDict()
    if counts:
        if not sfac_order:
            # 若未能读取到 SFAC 顺序，则按元素名称排序
            for elem in sorted(counts.keys()):
                ordered_counts[elem] = counts[elem]
        else:
            added = set()
            for elem in sfac_order:
                if elem in counts and elem not in added:
                    ordered_counts[elem] = counts[elem]
                    added.add(elem)
            # 补充在 SFAC 中未出现但统计到的元素
            for elem in counts:
                if elem not in ordered_counts:
                    ordered_counts[elem] = counts[elem]
    return ordered_counts


def composition_signature(ordered_counts: OrderedDict[str, int]) -> str:
    if not ordered_counts:
        return "未知配比"
    gcd_value = gcd_many(ordered_counts.values())
    parts: List[str] = []
    for elem, count in ordered_counts.items():
        base = count // gcd_value
        if base == 1:
            parts.append(elem)
        else:
            parts.append(f"{elem}{base}")
    return "".join(parts)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="统计指定迭代目录（IT*）中各个元素配比的数量。需在 PD 目录内运行。"
    )
    parser.add_argument(
        "iteration",
        help="迭代编号，例如 0 或 IT0。",
    )
    parser.add_argument(
        "--output",
        help="保存统计结果的文件路径（默认为直接打印）。",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pd_dir = Path.cwd()
    if pd_dir.name != "PD":
        print("警告：当前不在 PD 目录，建议切换到 PD 再运行。", file=sys.stderr)

    try:
        iteration = parse_iteration_label(args.iteration)
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    it_dir = pd_dir / iteration
    if not it_dir.is_dir():
        print(f"错误：未找到迭代目录 {it_dir}", file=sys.stderr)
        return 1

    res_files = sorted(it_dir.glob("*.res"))
    total_files = len(res_files)

    composition_counter: Counter[str] = Counter()

    show_progress = total_files > 0 and sys.stdout.isatty()
    bar_width = 40
    for idx, res_path in enumerate(res_files, start=1):
        ordered_counts = parse_composition(res_path)
        signature = composition_signature(ordered_counts)
        composition_counter[signature] += 1
        if show_progress:
            filled = int(bar_width * idx / total_files)
            bar = "#" * filled + "." * (bar_width - filled)
            sys.stdout.write(f"\r解析进度：[{bar}] {idx}/{total_files}")
            sys.stdout.flush()

    if show_progress:
        sys.stdout.write("\n")

    lines: List[str] = []
    lines.append(f"统计目录：{iteration}")
    lines.append(f"文件总数：{total_files}")
    lines.append("配比\t数量")
    for signature, count in composition_counter.most_common():
        lines.append(f"{signature}\t{count}")

    output_text = "\n".join(lines)
    if args.output:
        output_path = Path(args.output)
        if not output_path.is_absolute():
            output_path = (pd_dir / output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(output_text + "\n", encoding="utf-8")
        print(f"统计结果已写入：{output_path}")
    else:
        print(output_text)

    return 0


if __name__ == "__main__":
    sys.exit(main())
