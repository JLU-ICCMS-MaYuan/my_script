#!/usr/bin/env python3
"""
根据元素区间或区块(block)定义生成配比，并按分子式倍数扩展，可通过
--natom 限制最终原子总数。

示例：
    python generate_formula_ratios.py --atom Ce:1-10 Mg:1-10 H:1-30 -fu 1 2
    python generate_formula_ratios.py --block CeH2 MgH3 H --natom 80 -fu 1 2 3
    python generate_formula_ratios.py --atom Ce:1-3 Mg:1-5 H:1-10 --natom 40
"""

from __future__ import annotations

import argparse
import itertools
import re
import sys
from collections import Counter
from dataclasses import dataclass
from functools import reduce
from math import gcd
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

DEFAULT_OUTPUT = "composition.dat"
DEFAULT_NUM_STRUCTURES = 1


@dataclass
class GenerationConfig:
    """承载生成配比时需要的参数集合。"""

    mode: str
    symbols: List[str]
    display_symbols: List[str]
    display_indices: List[int]
    formula_units: Sequence[int]
    deduplicate: bool
    num_structures: int
    spec_map: Sequence[Tuple[str, range]]
    blocks: Sequence[Counter[str]]
    block_atom_limit: Optional[int]
    mh_range: Optional[Tuple[float, float]]
    max_atoms: Optional[int]
    reducible_divisors: Optional[Tuple[int, ...]]


def parse_element_spec(spec: str) -> Tuple[str, range]:
    """解析类似 Ce:1-3 的区间说明。"""

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
    """对计数做最大公因数约化，避免等比例重复。"""

    divisor = reduce(gcd, counts)
    divisor = divisor or 1
    return tuple(count // divisor for count in counts)


def format_composition(
    symbols: Sequence[str],
    counts: Sequence[int],
    display_symbols: Sequence[str],
    display_indices: Sequence[int],
) -> str:
    """根据显示顺序拼接化学式字符串。"""

    pairs = (
        (display_symbols[i], counts[display_indices[i]])
        for i in range(len(display_symbols))
    )
    parts = [f"{symbol}{count}" for symbol, count in pairs if count > 0]
    return "".join(parts)


def parse_block_spec(spec: str) -> Counter[str]:
    """解析区块表示（例如 CeH2）。"""

    pattern = re.compile(r"([A-Z][a-z]*)(\d*)")
    pos = 0
    counts: Counter[str] = Counter()
    spec = spec.strip()
    for match in pattern.finditer(spec):
        if match.start() != pos:
            raise ValueError(f"无法解析区块: '{spec}'")
        symbol, num = match.groups()
        value = int(num) if num else 1
        if value <= 0:
            raise ValueError(f"区块中包含非正整数: '{spec}'")
        counts[symbol] += value
        pos = match.end()
    if pos != len(spec) or not counts:
        raise ValueError(f"无法解析区块: '{spec}'")
    return counts


def iter_counts_from_atoms(spec_map: Sequence[Tuple[str, range]]) -> Iterator[Tuple[int, ...]]:
    """依据元素区间做笛卡尔积得到所有计数组合。"""

    ranges = [rng for _, rng in spec_map]
    for counts in itertools.product(*ranges):
        yield tuple(counts)


def iter_counts_from_blocks(
    blocks: Sequence[Counter[str]],
    symbols: Sequence[str],
    natom: int,
) -> Iterator[Tuple[int, ...]]:
    """在不超过 natom 的条件下按区块组合生成计数。"""

    index = {symbol: idx for idx, symbol in enumerate(symbols)}
    totals = [sum(counter.values()) for counter in blocks]
    counts = [0] * len(symbols)

    def backtrack(idx: int, remaining: int) -> Iterator[Tuple[int, ...]]:
        if idx == len(blocks):
            if any(counts):
                yield tuple(counts)
            return
        block_counts = blocks[idx]
        block_atoms = totals[idx]
        max_units = remaining // block_atoms if block_atoms else 0
        for qty in range(max_units + 1):
            if qty:
                for symbol, value in block_counts.items():
                    counts[index[symbol]] += value * qty
            yield from backtrack(idx + 1, remaining - block_atoms * qty)
            if qty:
                for symbol, value in block_counts.items():
                    counts[index[symbol]] -= value * qty

    yield from backtrack(0, natom)


def expand_formulas(
    symbols: Sequence[str],
    display_symbols: Sequence[str],
    display_indices: Sequence[int],
    base_counts: Iterable[Tuple[int, ...]],
    fus: Sequence[int],
    deduplicate: bool,
) -> Iterable[Tuple[str, Tuple[int, ...]]]:
    """对基础配比做约化、倍增与去重，输出最终化学式与计数。"""

    seen = set()
    for counts in base_counts:
        if not any(counts):
            continue
        reduced = reduce_counts(counts)
        if not any(reduced):
            continue
        for fu in fus:
            final_counts = tuple(value * fu for value in reduced)
            if deduplicate and final_counts in seen:
                continue
            seen.add(final_counts)
            formula = format_composition(symbols, final_counts, display_symbols, display_indices)
            if not formula:
                continue
            yield formula, final_counts


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="生成元素配比列表，可使用原子区间或区块(block)描述。",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--atom",
        nargs="+",
        metavar="SPEC",
        help="元素范围，格式如 Ce:1-3 Mg:1-3 H:1-20",
    )
    group.add_argument(
        "--block",
        nargs="+",
        metavar="BLOCK",
        help="区块描述，如 CeH2 MgH3 H",
    )
    parser.add_argument(
        "--natom",
        type=int,
        help="最大原子数限制；区块模式下表示基础配比所允许的原子总数",
    )
    parser.add_argument(
        "-fu",
        "--formula-unit",
        nargs="+",
        type=int,
        default=[1],
        help="分子式倍数，可写多个，如 -fu 1 2 3 表示生成1×,2×,3×分子式",
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
    parser.add_argument(
        "-f", "--mh-range",
        nargs=2,
        metavar=("LOW", "HIGH"),
        help="筛选金属:氢比例范围，例如 -f 1/6 1/8 表示金属总数/氢总数在 [1/6, 1/8] 之间",
    )
    parser.add_argument(
        "-od",
        "--order",
        nargs="+",
        metavar="ELEMENT",
        help="指定输出配比的元素顺序，例如 -od Ce Mg H",
    )
    parser.add_argument(
        "--require-reducible",
        nargs="+",
        type=int,
        metavar="DIVISOR",
        help="仅保留所有元素计数可同时被指定因子整除的配比，例如 --require-reducible 4 6 8",
    )
    return parser


def parse_mh_range(parser: argparse.ArgumentParser, raw_range: Optional[Sequence[str]]) -> Optional[Tuple[float, float]]:
    """解析金属:氢比例上下限。"""

    if not raw_range:
        return None

    def parse_ratio(text: str) -> float:
        if "/" in text:
            numerator, denominator = text.split("/", 1)
            return float(numerator) / float(denominator)
        return float(text)

    try:
        low = parse_ratio(raw_range[0])
        high = parse_ratio(raw_range[1])
    except Exception:
        parser.error("比例范围参数 -f 必须是形如 1/6 或 0.2 的数字")
    if low > high:
        parser.error("-f/--mh-range 的下限不能大于上限")
    return (low, high)


def build_config(parser: argparse.ArgumentParser, args: argparse.Namespace) -> GenerationConfig:
    """校验参数并打包成 GenerationConfig。"""

    if any(fu <= 0 for fu in args.formula_unit):
        parser.error("所有 -fu 参数必须为正整数")
    if args.num_structures <= 0:
        parser.error("每个配比的结构数 (-n/--num-structures) 必须为正整数")
    if args.natom is not None and args.natom <= 0:
        parser.error("--natom 需要是正整数")
    reducible_divisors: Optional[Tuple[int, ...]] = None
    if args.require_reducible:
        if any(divisor <= 1 for divisor in args.require_reducible):
            parser.error("--require-reducible 需要大于 1 的正整数")
        reducible_divisors = tuple(args.require_reducible)

    mode: Optional[str] = None
    spec_map: List[Tuple[str, range]] = []
    block_specs: List[Counter[str]] = []
    symbols: List[str] = []

    if args.atom:
        try:
            spec_map = [parse_element_spec(spec) for spec in args.atom]
        except ValueError as exc:
            parser.error(str(exc))
        symbols = [symbol for symbol, _ in spec_map]
        mode = "atom"
    elif args.block:
        if args.natom is None:
            parser.error("区块模式下必须设置 --natom 作为原子数上限")
        seen_symbols: List[str] = []
        try:
            for block in args.block:
                counts = parse_block_spec(block)
                block_specs.append(counts)
                for sym in counts:
                    if sym not in seen_symbols:
                        seen_symbols.append(sym)
        except ValueError as exc:
            parser.error(str(exc))
        symbols = seen_symbols
        mode = "block"
    else:  # pragma: no cover
        parser.error("必须指定 --atom 或 --block")

    if args.order:
        order = args.order
        if len(order) != len(symbols):
            parser.error("-od/--order 需要指定与元素数量相同的符号")
        if sorted(order) != sorted(symbols):
            parser.error("-od/--order 必须包含所有元素且仅包含这些元素")
        display_symbols = order
    else:
        display_symbols = symbols
    display_indices = [symbols.index(sym) for sym in display_symbols]

    mh_range = parse_mh_range(parser, args.mh_range)
    max_atoms = args.natom if args.natom is not None else None

    return GenerationConfig(
        mode=mode,
        symbols=symbols,
        display_symbols=display_symbols,
        display_indices=display_indices,
        formula_units=args.formula_unit,
        deduplicate=not args.allow_duplicates,
        num_structures=args.num_structures,
        spec_map=spec_map,
        blocks=block_specs,
        block_atom_limit=args.natom,
        mh_range=mh_range,
        max_atoms=max_atoms,
        reducible_divisors=reducible_divisors,
    )


def build_base_iterator(config: GenerationConfig) -> Iterable[Tuple[int, ...]]:
    """根据模式返回基础配比迭代器。"""

    if config.mode == "atom":
        return iter_counts_from_atoms(config.spec_map)
    if config.block_atom_limit is None:
        raise ValueError("区块模式缺少 natom 限制")
    return iter_counts_from_blocks(config.blocks, config.symbols, config.block_atom_limit)


def counts_in_mh_range(config: GenerationConfig, counts: Sequence[int]) -> bool:
    """判断金属:氢比例是否满足限制。"""

    if not config.mh_range:
        return True
    metal = sum(c for s, c in zip(config.symbols, counts) if s != "H")
    hydrogen = sum(c for s, c in zip(config.symbols, counts) if s == "H")
    if hydrogen == 0:
        return False
    ratio = metal / hydrogen
    low, high = config.mh_range
    return low <= ratio <= high


def counts_reducible(counts: Sequence[int], divisors: Optional[Tuple[int, ...]]) -> bool:
    """判断计数是否可被任一因子整体约分。"""

    if not divisors:
        return True
    return any(all(count % divisor == 0 for count in counts) for divisor in divisors)


def filter_formulas(
    config: GenerationConfig,
    base_iter: Iterable[Tuple[int, ...]],
) -> Iterator[str]:
    """结合所有限制条件，输出最终的化学式字符串。"""

    for formula, final_counts in expand_formulas(
        config.symbols,
        config.display_symbols,
        config.display_indices,
        base_iter,
        config.formula_units,
        deduplicate=config.deduplicate,
    ):
        if config.mode == "block" and any(count == 0 for count in final_counts):
            continue
        if config.max_atoms is not None and sum(final_counts) > config.max_atoms:
            continue
        if not counts_in_mh_range(config, final_counts):
            continue
        if not counts_reducible(final_counts, config.reducible_divisors):
            continue
        yield formula


def write_output(path: str, formulas: Sequence[str], num_structures: int) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for formula in formulas:
            handle.write(f"{formula} {num_structures}\n")


def main(argv: Sequence[str]) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    config = build_config(parser, args)

    base_iter = build_base_iterator(config)
    formulas = list(filter_formulas(config, base_iter))

    write_output(args.output, formulas, config.num_structures)
    print(f"# 已写入 {args.output}，共 {len(formulas)} 条配比。")
    if config.max_atoms is not None and not formulas:
        print("提示：当前最大原子数限制可能过小，可尝试调大 --natom。")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
