#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
import itertools
import shlex
import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from functools import reduce
from math import gcd
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

# --- Data Section for USPEX-like method ---
# This dictionary contains empirical atomic volumes at 0 GPa (in Angstrom^3/atom) for the first 86 elements.
EMPIRICAL_ATOMIC_VOLUMES = {
    # Period 1
    'H': 5.8, 'He': 20.0,
    # Period 2
    'Li': 21.6, 'Be': 8.1, 'B': 7.5, 'C': 5.7, 'N': 9.2, 'O': 6.9, 'F': 5.5, 'Ne': 22.1,
    # Period 3
    'Na': 39.4, 'Mg': 23.2, 'Al': 16.6, 'Si': 20.0, 'P': 19.5, 'S': 23.9, 'Cl': 20.5, 'Ar': 37.0,
    # Period 4
    'K': 75.3, 'Ca': 43.2, 'Sc': 25.0, 'Ti': 17.7, 'V': 13.8, 'Cr': 12.0, 'Mn': 12.1,
    'Fe': 11.7, 'Co': 11.1, 'Ni': 10.9, 'Cu': 11.8, 'Zn': 15.2, 'Ga': 19.5, 'Ge': 22.6,
    'As': 21.4, 'Se': 27.2, 'Br': 29.0, 'Kr': 48.2,
    # Period 5
    'Rb': 92.2, 'Sr': 56.4, 'Y': 33.0, 'Zr': 23.2, 'Nb': 17.9, 'Mo': 15.6, 'Tc': 14.1,
    'Ru': 13.5, 'Rh': 13.7, 'Pd': 14.7, 'Ag': 17.0, 'Cd': 21.5, 'In': 26.1, 'Sn': 27.0,
    'Sb': 30.2, 'Te': 33.8, 'I': 36.9, 'Xe': 59.8,
    # Period 6
    'Cs': 116.3, 'Ba': 63.9,
    # Lanthanides
    'La': 37.3, 'Ce': 34.5, 'Pr': 34.4, 'Nd': 34.2, 'Pm': 33.9, 'Sm': 32.9, 'Eu': 48.0,
    'Gd': 32.4, 'Tb': 31.6, 'Dy': 31.2, 'Ho': 30.8, 'Er': 30.5, 'Tm': 30.1, 'Yb': 41.3, 'Lu': 29.5,
    # Rest of Period 6
    'Hf': 22.2, 'Ta': 18.0, 'W': 15.8, 'Re': 14.7, 'Os': 14.0, 'Ir': 14.1, 'Pt': 15.1,
    'Au': 16.9, 'Hg': 23.0, 'Tl': 28.5, 'Pb': 30.3, 'Bi': 35.3, 'Po': 35.8, 'At': 39.0, 'Rn': 50.5
}

@dataclass
class GenerationConfig:
    mode: str
    symbols: List[str]
    display_symbols: List[str]
    display_indices: List[int]
    formula_units: Sequence[int]
    deduplicate: bool
    num_structures: int
    spec_map: Sequence[Tuple[str, range]]
    blocks: Sequence[Counter]
    block_atom_limit: Optional[int]
    mh_range: Optional[Tuple[float, float]]
    max_atoms: Optional[int]
    reducible_divisors: Optional[Tuple[int, ...]]


def setup_arg_parser():
    """Sets up the argument parser for the script with short aliases."""
    parser = argparse.ArgumentParser(
        description="Generate CALYPSO input.dat files for a range of compositions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # --- Core Arguments ---
    parser.add_argument(
        '-e', '--elements', nargs='+', required=True, type=str,
        help="List of element names.\nExample: -e Si O"
    )
    parser.add_argument(
        '-r', '--radii', nargs='+', required=True, type=float,
        help="List of covalent radii for each element (in Bohr).\nRequired for distance matrix.\nExample: -r 1.11 0.60"
    )
    parser.add_argument(
        '-p', '--popsize', required=True, type=int,
        help="Population size for the CALYPSO run (e.g., 30)."
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
        "-n",
        "--num-structures",
        type=int,
        default=1,
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
    # --- Volume Calculation Method Arguments ---
    parser.add_argument(
        '-m', '--method', type=str, choices=['radius', 'uspex'], default='radius',
        help="Method for volume calculation: 'radius' (default) or 'uspex'."
    )
    parser.add_argument(
        '-P', '--pressure', type=float,
        help="Target pressure in GPa (required for '--method uspex').\n(Note: Uppercase 'P' to avoid conflict)"
    )
    return parser


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


def parse_block_spec(spec: str) -> Counter:
    """解析区块表示（例如 CeH2）。"""
    pattern = re.compile(r"([A-Z][a-z]*)(\d*)")
    pos = 0
    counts: Counter = Counter()
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
    blocks: Sequence[Counter],
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


def format_composition(display_symbols: Sequence[str], display_indices: Sequence[int], counts: Sequence[int]) -> str:
    pairs = (
        (display_symbols[i], counts[display_indices[i]])
        for i in range(len(display_symbols))
    )
    parts = [f"{symbol}{count}" for symbol, count in pairs if count > 0]
    return "".join(parts)


def expand_formulas(
    base_counts: Iterable[Tuple[int, ...]],
    fus: Sequence[int],
    deduplicate: bool,
) -> Iterable[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """对基础配比做约化、倍增与去重，输出最终计数与约化计数。"""
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
            yield final_counts, reduced


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


def counts_in_mh_range(symbols: Sequence[str], counts: Sequence[int], mh_range: Optional[Tuple[float, float]]) -> bool:
    """判断金属:氢比例是否满足限制。"""
    if not mh_range:
        return True
    metal = sum(c for s, c in zip(symbols, counts) if s != "H")
    hydrogen = sum(c for s, c in zip(symbols, counts) if s == "H")
    if hydrogen == 0:
        return False
    ratio = metal / hydrogen
    low, high = mh_range
    return low <= ratio <= high


def counts_reducible(counts: Sequence[int], divisors: Optional[Tuple[int, ...]]) -> bool:
    """判断计数是否可被任一因子整体约分。"""
    if not divisors:
        return True
    return any(all(count % divisor == 0 for count in counts) for divisor in divisors)


def build_config(parser: argparse.ArgumentParser, args: argparse.Namespace) -> GenerationConfig:
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
    block_specs: List[Counter] = []
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

def calculate_volume_radius(radii_bohr, composition):
    """Estimates the cell volume based on covalent radii (original method)."""
    scaled_radii_angstrom = radii_bohr * 0.529177 * 0.7
    sphere_volumes = (4/3) * np.pi * np.power(scaled_radii_angstrom, 3)
    total_volume = np.sum(sphere_volumes * composition)
    return total_volume * 1.1

def calculate_volume_uspex(elements, composition, pressure_gpa):
    """Estimates the cell volume using the USPEX-like empirical method."""
    v0 = 0.0
    for el, num in zip(elements, composition):
        if el not in EMPIRICAL_ATOMIC_VOLUMES:
            raise ValueError(f"Error: Empirical volume for element '{el}' is not in the database.")
        v0 += EMPIRICAL_ATOMIC_VOLUMES[el] * num
    
    if pressure_gpa is None or pressure_gpa == 0:
        return v0
    
    a = 0.08  # Empirical constant for pressure response
    b = 2.5   # Empirical constant for curvature
    volume_p = v0 * (1.0 / (1.0 + a * pressure_gpa))**(1.0/b)
    return volume_p

def calculate_distance_matrix(radii_bohr):
    """Calculates the NxN matrix of interatomic distances from radii."""
    distance_matrix_angstrom = (radii_bohr[:, np.newaxis] + radii_bohr) * 0.529177 * 0.7
    return distance_matrix_angstrom

def generate_input_file(elements, radii, composition, popsize, num_elements, method, pressure):
    """Generates a single input.dat file for a given composition."""
    # --- Volume Calculation ---
    if method == 'radius':
        volume = calculate_volume_radius(radii, composition)
    elif method == 'uspex':
        volume = calculate_volume_uspex(elements, composition, pressure)
    
    # --- Distance Calculation ---
    distance_matrix = calculate_distance_matrix(radii)

    # --- Directory and File Setup ---
    dir_name = "".join([f"{el}{num}" for el, num in zip(elements, composition)])
    os.makedirs(dir_name, exist_ok=True)
    output_path = os.path.join(dir_name, 'input.dat')
    
    # --- File Content Generation ---
    system_name = "-".join(elements)
    atom_names = " ".join(elements)
    atom_counts = " ".join(map(str, composition))
    distance_lines = [" ".join([f"{dist:.2f}" for dist in row]) for row in distance_matrix]

    # --- Writing to input.dat ---
    with open(output_path, 'w') as f:
        f.write(f"SystemName = {system_name}\n")
        f.write(f"NumberOfSpecies = {num_elements}\n")
        f.write(f"NameOfAtoms = {atom_names}\n")
        f.write(f"NumberOfAtoms = {atom_counts}\n")
        f.write("NumberOfFormula = 1 1\n")
        f.write(f"Volume = {volume:.4f}\n\n")
        f.write("@DistanceOfIon\n")
        f.write("\n".join(distance_lines) + "\n")
        f.write("@End\n\n")
        f.write(f"PopSize = {popsize}\n")
        f.write("SpeSpaceGroup = 8 230\n")
        f.write("Split = T\n")
        f.write("VSC = F\n")
        f.write("MaxNumAtom = 80\n")
        f.write("Command = sh submit.sh\n")
        f.write("MaxStep = 10\n")


def build_base_iterator(config: GenerationConfig) -> Iterable[Tuple[int, ...]]:
    """根据模式返回基础配比迭代器。"""
    if config.mode == "atom":
        return iter_counts_from_atoms(config.spec_map)
    if config.block_atom_limit is None:
        raise ValueError("区块模式缺少 natom 限制")
    return iter_counts_from_blocks(config.blocks, config.symbols, config.block_atom_limit)


def filter_compositions(config: GenerationConfig) -> Iterable[Tuple[int, ...]]:
    """结合所有限制条件，输出最终计数组合。"""
    base_iter = build_base_iterator(config)
    for final_counts, _ in expand_formulas(
        base_iter,
        config.formula_units,
        deduplicate=config.deduplicate,
    ):
        if config.mode == "block" and any(count == 0 for count in final_counts):
            continue
        if config.max_atoms is not None and sum(final_counts) > config.max_atoms:
            continue
        if not counts_in_mh_range(config.symbols, final_counts, config.mh_range):
            continue
        if not counts_reducible(final_counts, config.reducible_divisors):
            continue
        yield final_counts


def format_formula_for_log(config: GenerationConfig, counts: Sequence[int]) -> str:
    return format_composition(config.display_symbols, config.display_indices, counts)

def write_parameter_log(args, filename="paras.dat"):
    """Persist the latest CLI arguments for easier reruns."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        need_separator = os.path.exists(filename) and os.path.getsize(filename) > 0
        with open(filename, 'a') as f:
            if need_separator:
                f.write("\n")
            quoted_args = " ".join(shlex.quote(arg) for arg in sys.argv)
            f.write(f"[{timestamp}] Invocation: {quoted_args}\n")
            for key, value in sorted(vars(args).items()):
                if isinstance(value, (list, tuple)):
                    value_str = " ".join(map(str, value))
                else:
                    value_str = str(value)
                f.write(f"{key} = {value_str}\n")
            f.write("----\n")
    except OSError as e:
        print(f"Warning: 无法写入参数记录文件 '{filename}': {e}")

def main():
    """Main function to parse arguments and generate files."""
    parser = setup_arg_parser()
    args = parser.parse_args()
    config = build_config(parser, args)

    # --- Argument Validation ---
    if args.method == 'uspex' and args.pressure is None:
        parser.error("--pressure (-P) is required when using --method uspex")

    num_elements = len(args.elements)
    if not (len(args.radii) == num_elements and len(config.symbols) == num_elements):
        print("Error: The number of arguments for --elements (-e), --radii (-r), and --atom/--block must be the same.")
        sys.exit(1)
    if args.elements != config.symbols:
        print("Error: --elements 的顺序必须与 --atom/--block 中元素顺序一致。")
        sys.exit(1)

    write_parameter_log(args)

    compositions = list(filter_compositions(config))
    total_files = len(compositions)
    print(f"Found {total_files} composition(s). Generating files using '{args.method}' method...")

    all_dir_names = []
    for i, comp in enumerate(compositions):
        current_composition = np.array(comp)
        dir_name = "".join([f\"{el}{num}\" for el, num in zip(args.elements, current_composition)])
        all_dir_names.append(dir_name)

        log_formula = format_formula_for_log(config, current_composition)
        log_text = log_formula if log_formula else dir_name
        print(f\"[{i+1}/{total_files}] Generating input for {log_text}...\")
        generate_input_file(
            elements=args.elements,
            radii=np.array(args.radii),
            composition=current_composition,
            popsize=args.popsize,
            num_elements=num_elements,
            method=args.method,
            pressure=args.pressure
        )

    try:
        output_filename = "composition.dat"
        with open(output_filename, 'w') as f:
            for name in all_dir_names:
                f.write(f\"./{name}\\n\")
        print(f\"\\nSuccessfully created all {total_files} input files.\")
        print(f\"All {len(all_dir_names)} directory paths have been written to '{output_filename}'.\")
    except IOError as e:
        print(f\"\\nWarning: Could not write to file '{output_filename}': {e}\")
        print(\"Input files were still created successfully.\")


if __name__ == "__main__":
    main()
