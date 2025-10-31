#!/usr/bin/env python3
"""
根据指定元素及其原子数范围生成配比，并化为最简形式。

示例：
    python generate_formula_ratios.py Ce:1-10 Mg:1-10 H:1-30 -fu 2
"""

from __future__ import annotations

import argparse
import itertools
import sys
from functools import reduce
from math import gcd
from typing import Iterable, List, Sequence, Tuple

DEFAULT_OUTPUT = "composition.dat"
DEFAULT_NUM_STRUCTURES = 1


def parse_element_spec(spec: str) -> Tuple[str, range]:
    """将形如 Ce:1-10 的参数解析为 (元素符号, range)。"""
    if ":" not in spec:
        raise ValueError(f"元素参数缺少冒号: '{spec}'")
    symbol, span = spec.split(":", 1)
    symbol = symbol.strip()
    if not symbol:
        raise ValueError(f"元素符号为空: '{spec}'")

    if "-" not in span:
        raise ValueError(f"元素范围缺少连字符: '{spec}'")
    start_str, end_str = span.split("-", 1)
    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError as exc:
        raise ValueError(f"范围必须是整数: '{spec}'") from exc

    if start <= 0 or end <= 0:
        raise ValueError(f"范围必须为正整数: '{spec}'")
    if start > end:
        raise ValueError(f"下限不能大于上限: '{spec}'")

    return symbol, range(start, end + 1)


def reduce_counts(counts: Sequence[int]) -> Tuple[int, ...]:
    """将计数约化为最简整数比。"""
    divisor = reduce(gcd, counts)
    return tuple(count // divisor for count in counts)


def format_composition(symbols: Sequence[str], counts: Sequence[int]) -> str:
    """按元素顺序输出紧凑格式形如 Ce1Mg2H3 的字符串。"""
    return "".join(f"{symbol}{count}" for symbol, count in zip(symbols, counts))


def generate_ratios(
    spec_map: Sequence[Tuple[str, range]],
    fu: int,
    deduplicate: bool = True,
) -> Iterable[Tuple[str, Tuple[int, ...]]]:
    """生成(文本行, 约化后数量) 对，以便后续扩展。"""
    symbols = [symbol for symbol, _ in spec_map]
    ranges = [r for _, r in spec_map]
    seen = set()

    for counts in itertools.product(*ranges):
        reduced = reduce_counts(counts)
        final_counts = tuple(count * fu for count in reduced)
        if deduplicate:
            key = final_counts
            if key in seen:
                continue
            seen.add(key)
        yield format_composition(symbols, final_counts), final_counts


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="生成元素配比列表，并可指定分子式倍数。",
    )
    parser.add_argument(
        "element_specs",
        nargs="+",
        help="元素范围，格式如 Ce:1-10 Mg:1-10 H:1-30",
    )
    parser.add_argument(
        "-fu",
        "--formula-unit",
        type=int,
        default=1,
        help="分子式倍数，默认 1 表示最简配比",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"输出文件路径（默认: {DEFAULT_OUTPUT}）",
    )
    parser.add_argument(
        "-n",
        "--num-structures",
        type=int,
        default=DEFAULT_NUM_STRUCTURES,
        help="每个配比对应的结构数量，默认 1",
    )
    parser.add_argument(
        "--allow-duplicates",
        action="store_true",
        help="保留重复配比（默认去重）",
    )
    return parser


def main(argv: Sequence[str]) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    if args.formula_unit <= 0:
        parser.error("公式倍数 (--formula-unit) 必须为正整数")

    try:
        spec_map = [parse_element_spec(spec) for spec in args.element_specs]
    except ValueError as exc:
        parser.error(str(exc))

    if args.num_structures <= 0:
        parser.error("每个配比的结构数 (-n/--num-structures) 必须为正整数")

    lines: List[str] = []
    lines.append(
        f"# 元素范围: {' '.join(args.element_specs)}; fu={args.formula_unit}; 去重={'否' if args.allow_duplicates else '是'}; 每配比结构数={args.num_structures}"
    )

    for line, _ in generate_ratios(
        spec_map,
        args.formula_unit,
        deduplicate=not args.allow_duplicates,
    ):
        lines.append(f"{line} {args.num_structures}")

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")

    print(f"# 已写入 {args.output}，共 {len(lines) - 1} 条配比。")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
