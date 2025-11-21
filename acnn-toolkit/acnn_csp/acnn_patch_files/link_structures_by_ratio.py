#!/usr/bin/env python3
"""按配比条件建立软链接的工具。

示例：
    python link_structures_by_ratio.py IT0 --ratio Ce:H<1:2 Mg:H<1:2
    python link_structures_by_ratio.py --ratio Ce:H>2:1 --dest selected --dry-run

该脚本需在 PD 目录运行，遍历指定 IT* 目录下的 .res 文件，
对约分后的化学式判断比例条件，满足全部条件的结构会在目标目录
创建指向原文件的软链接，方便进一步处理。
"""

from __future__ import annotations

import argparse
import math
import os
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

from ase.io import read  # type: ignore


RatioCond = Tuple[str, str, str, int, int]


@dataclass
class StructureRecord:
    iteration: str
    path: Path
    reduced_counts: Counter[str]


def parse_iteration_label(token: str) -> str:
    token = token.strip().upper()
    if not token:
        raise ValueError("迭代编号不能为空")
    if token.startswith("IT"):
        return token
    return f"IT{token}"


def discover_it_dirs(root: Path, iterations: Sequence[str]) -> List[Path]:
    if iterations:
        dirs: List[Path] = []
        for label in iterations:
            path = root / parse_iteration_label(label)
            if path.is_dir():
                dirs.append(path)
        return dirs
    return sorted(p for p in root.glob("IT*") if p.is_dir())


def reduce_counts(counts: Counter[str]) -> Counter[str]:
    gcd_val = 0
    for value in counts.values():
        gcd_val = math.gcd(gcd_val, value)
    gcd_val = max(gcd_val, 1)
    reduced = Counter()
    for elem, value in counts.items():
        reduced[elem] = value // gcd_val
    return reduced


def _parse_ratio_value(text: str) -> Fraction:
    parts = [seg.strip() for seg in text.split(":") if seg.strip()]
    if not parts:
        raise ValueError(f"比例值不能为空：{text}")
    try:
        if len(parts) == 1:
            ratio = Fraction(parts[0])
        elif len(parts) == 2:
            numerator = Fraction(parts[0])
            denominator = Fraction(parts[1])
            if denominator == 0:
                raise ValueError
            ratio = numerator / denominator
        else:
            raise ValueError
    except ValueError as exc:
        raise ValueError(f"比例值必须为数字：{text}") from exc
    if ratio <= 0:
        raise ValueError(f"比例值必须大于 0：{text}")
    return ratio


def _parse_element_pair(text: str) -> Tuple[str, str]:
    tokens = [seg.strip() for seg in text.split(":") if seg.strip()]
    if len(tokens) != 2:
        raise ValueError(f"元素组合格式错误：{text}")
    return tokens[0], tokens[1]


def parse_ratio_token(token: str) -> List[RatioCond]:
    comparators = ("<=", ">=", "<", ">", "=")
    limits: List[RatioCond] = []
    for part in (seg.strip() for seg in token.split(",")):
        if not part:
            continue
        matched = False
        for comparator in comparators:
            if comparator in part:
                left, right = part.split(comparator, 1)
                elem1, elem2 = _parse_element_pair(left)
                ratio = _parse_ratio_value(right)
                limits.append((elem1, elem2, comparator, ratio.numerator, ratio.denominator))
                matched = True
                break
        if matched:
            continue
        segments = part.split(":")
        if len(segments) != 3:
            raise ValueError(f"比例格式错误：{part}")
        elem1, elem2, ratio_part = (s.strip() for s in segments)
        ratio = _parse_ratio_value(ratio_part)
        limits.append((elem1, elem2, "<=", ratio.numerator, ratio.denominator))
    if not limits:
        raise ValueError(f"比例格式错误：{token}")
    return limits


def match_ratios(record: StructureRecord, ratio_conditions: Sequence[RatioCond]) -> bool:
    if not ratio_conditions:
        return False
    reduced = record.reduced_counts
    for elem1, elem2, comparator, num, den in ratio_conditions:
        v1 = reduced.get(elem1, 0)
        v2 = reduced.get(elem2, 0)
        if v1 == 0 or v2 == 0:
            return False
        left = v1 * den
        right = v2 * num
        if comparator == "<" and not (left < right):
            return False
        if comparator == "<=" and not (left <= right):
            return False
        if comparator == ">" and not (left > right):
            return False
        if comparator == ">=" and not (left >= right):
            return False
        if comparator == "=" and not (left == right):
            return False
    return True


def load_records(it_dirs: Iterable[Path]) -> List[StructureRecord]:
    records: List[StructureRecord] = []
    for it_dir in it_dirs:
        for res_path in sorted(it_dir.glob("*.res")):
            try:
                atoms = read(res_path, format="res")
            except Exception as exc:  # pragma: no cover - 解析失败
                print(f"警告：{res_path} 解析失败 ({exc})")
                continue
            counts = Counter(atoms.get_chemical_symbols())
            reduced = reduce_counts(counts)
            records.append(
                StructureRecord(iteration=it_dir.name, path=res_path, reduced_counts=reduced)
            )
    return records


def create_symlink(src: Path, dst: Path, dry_run: bool) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if dry_run:
        print(f"[模拟] {dst} -> {src}")
        return
    dst.symlink_to(src.resolve())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="按配比条件建立软链接")
    parser.add_argument(
        "iterations",
        nargs="*",
        help="待扫描的迭代编号，例如 IT0 或 0；缺省时扫描所有 IT*。",
    )
    parser.add_argument(
        "--ratio",
        nargs="+",
        action="extend",
        metavar="COND",
        required=True,
        help="比例条件，形如 Ce:H<1:2；一次可提供多个条件，视为 AND 关系。",
    )
    parser.add_argument(
        "--dest",
        default="linked_structs",
        help="软链接输出目录，默认 linked_structs",
    )
    parser.add_argument("--dry-run", action="store_true", help="仅打印匹配结果，不实际建立链接")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    root = Path.cwd()
    if root.name != "PD":
        print("⚠️  建议在 PD 目录运行，当前目录：", root)

    try:
        it_dirs = discover_it_dirs(root, args.iterations)
    except ValueError as exc:
        parser.error(str(exc))

    if not it_dirs:
        parser.error("未找到任何 IT* 目录")

    ratio_conditions: List[RatioCond] = []
    try:
        for token in args.ratio:
            ratio_conditions.extend(parse_ratio_token(token))
    except ValueError as exc:
        parser.error(str(exc))

    records = load_records(it_dirs)
    if not records:
        parser.error("未解析到任何结构")

    matched = [rec for rec in records if match_ratios(rec, ratio_conditions)]
    dest_root = root / args.dest

    print(f"共解析 {len(records)} 个结构，命中 {len(matched)} 个。")
    if not matched:
        return 0

    dest_root.mkdir(parents=True, exist_ok=True)
    for rec in matched:
        link_name = f"{rec.iteration}_{rec.path.name}"
        link_path = dest_root / link_name
        create_symlink(rec.path, link_path, args.dry_run)

    if args.dry_run:
        print("已在模拟模式下列出所有软链接。")
    else:
        print(f"软链接已创建于 {dest_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
