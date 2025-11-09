#!/usr/bin/env python3
"""
基于 ASE 的配比筛选工具。

核心功能：
1. 读取指定 IT* 目录下的 .res 结构文件，统计配比分布。
2. 按元素比例约束或指定配比签名筛选 / 剔除结构。
3. 支持将命中的结构移动到垃圾桶目录，或恢复垃圾桶中的结构。

使用示例：
    python filter_composition_for_PD.py IT0 --stats
    python filter_composition_for_PD.py --ratio-limit Ce:Mg:5 --filter
    python filter_composition_for_PD.py --ratio-limit Ce:Mg:5 -c Ce9Mg8H1 --move-invalid
    python filter_composition_for_PD.py --restore --yes
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from ase.io import read  # type: ignore

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    def tqdm(iterable, *args, **kwargs):
        return iterable


@dataclass
class StructureRecord:
    path: Path
    iteration: str
    composition: Counter[str]
    signature: str
    canonical_key: Tuple[Tuple[str, int], ...]
    unique_elements: int


# -----------------------------------------------------------------------------
# 基础工具
# -----------------------------------------------------------------------------


def parse_iteration_label(label: str) -> str:
    token = label.strip()
    if not token:
        raise ValueError("迭代编号不能为空")
    upper = token.upper()
    if upper.startswith("IT"):
        return upper
    return f"IT{upper}"


def discover_it_dirs(root: Path, iteration: Optional[str]) -> List[Path]:
    if iteration:
        target = root / iteration
        if not target.is_dir():
            return []
        return [target]
    return sorted(p for p in root.glob("IT*") if p.is_dir())


def available_workers() -> int:
    cpu = os.cpu_count() or 1
    return max(1, cpu - 1)


def reduce_counts(counts: Counter[str]) -> Counter[str]:
    gcd_val = 0
    for value in counts.values():
        gcd_val = math.gcd(gcd_val, value) if gcd_val else value
    gcd_val = max(gcd_val, 1)
    reduced = Counter()
    for elem, value in counts.items():
        reduced[elem] = value // gcd_val
    return reduced


def counts_to_signature(counts: Counter[str]) -> str:
    if not counts:
        return "未知配比"
    parts = []
    for elem in sorted(counts):
        value = counts[elem]
        parts.append(f"{elem}{value if value > 1 else ''}")
    return "".join(parts)


def counts_to_key(counts: Counter[str]) -> Tuple[Tuple[str, int], ...]:
    reduced = reduce_counts(counts)
    return tuple(sorted((elem, reduced[elem]) for elem in reduced))


def parse_signature_to_key(signature: str) -> Tuple[Tuple[str, int], ...]:
    if not signature.strip():
        return tuple()
    tokens: Counter[str] = Counter()
    elem = ""
    num = ""
    for ch in signature.strip():
        if ch.isalpha():
            if ch.isupper():
                if elem:
                    value = int(num) if num else 1
                    tokens[elem] += value
                    num = ""
                elem = ch
            else:
                if elem:
                    elem += ch
        elif ch.isdigit():
            num += ch
        else:
            continue
    if elem:
        value = int(num) if num else 1
        tokens[elem] += value
    return counts_to_key(tokens)


def _read_res_file(path_str: str) -> Tuple[Optional[Counter[str]], Optional[str]]:
    path = Path(path_str)
    try:
        atoms = read(path, format="res")
    except Exception as exc:  # pragma: no cover - ASE 解析异常
        return None, f"{path}: 解析失败 ({exc})"
    counts = Counter(atoms.get_chemical_symbols())
    return counts, None


def _move_file_task(src: Path, dst: Path) -> Optional[str]:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
        return None
    except Exception as exc:  # pragma: no cover - IO 错误
        return str(exc)


def parse_ratio_limit(token: str) -> Tuple[str, str, float]:
    parts = token.split(":")
    if len(parts) != 3:
        raise ValueError(f"比例约束格式错误：{token}")
    elem1, elem2, ratio = (p.strip() for p in parts)
    if not elem1 or not elem2:
        raise ValueError(f"元素名称不能为空：{token}")
    try:
        value = float(ratio)
    except ValueError as exc:
        raise ValueError(f"比例值必须为数字：{token}") from exc
    if value <= 0:
        raise ValueError(f"比例值必须大于 0：{token}")
    return elem1, elem2, value


def format_ratio(counts: Counter[str], elem1: str, elem2: str) -> str:
    n1 = counts.get(elem1, 0)
    n2 = counts.get(elem2, 0)
    if n1 == 0 or n2 == 0:
        return f"{elem1}:{elem2} = N/A"
    gcd_val = math.gcd(n1, n2)
    return f"{elem1}:{elem2} = {n1 // gcd_val}:{n2 // gcd_val}"


# -----------------------------------------------------------------------------
# 核心逻辑
# -----------------------------------------------------------------------------


def load_structures(it_dirs: Iterable[Path]) -> Tuple[List[StructureRecord], List[str]]:
    entries: List[Tuple[str, Path]] = []
    for it_dir in it_dirs:
        for res_path in sorted(it_dir.glob("*.res")):
            entries.append((it_dir.name, res_path))

    if not entries:
        return [], []

    records: List[StructureRecord] = []
    warnings: List[str] = []
    workers = available_workers()
    use_parallel = len(entries) > 1 and workers > 1

    if use_parallel:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_read_res_file, str(path)): (iteration, path)
                for iteration, path in entries
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="解析结构", unit="个"):
                iteration, res_path = futures[future]
                counts, error = future.result()
                if error or counts is None:
                    if error:
                        warnings.append(error)
                    continue
                records.append(
                    StructureRecord(
                        path=res_path,
                        iteration=iteration,
                        composition=counts,
                        signature=counts_to_signature(reduce_counts(counts)),
                        canonical_key=counts_to_key(counts),
                        unique_elements=len(counts),
                    )
                )
    else:
        for iteration, res_path in tqdm(entries, desc="解析结构", unit="个"):
            counts, error = _read_res_file(str(res_path))
            if error or counts is None:
                if error:
                    warnings.append(error)
                continue
            records.append(
                StructureRecord(
                    path=res_path,
                    iteration=iteration,
                    composition=counts,
                    signature=counts_to_signature(reduce_counts(counts)),
                    canonical_key=counts_to_key(counts),
                    unique_elements=len(counts),
                )
            )

    return records, warnings


def should_skip_by_element_count(record: StructureRecord, exclude_counts: Sequence[int]) -> bool:
    return bool(exclude_counts) and record.unique_elements in exclude_counts


def check_ratio_constraints(
    record: StructureRecord,
    ratio_limits: Sequence[Tuple[str, str, float]],
) -> Tuple[bool, List[str]]:
    if not ratio_limits:
        return True, []
    violations: List[str] = []
    for elem1, elem2, max_ratio in ratio_limits:
        n1 = record.composition.get(elem1, 0)
        n2 = record.composition.get(elem2, 0)
        if n1 == 0 or n2 == 0:
            violations.append(f"{elem1} 或 {elem2} 缺失")
            continue
        ratio = max(n1, n2) / min(n1, n2)
        if ratio > max_ratio:
            violations.append(f"{format_ratio(record.composition, elem1, elem2)} > {max_ratio}:1")
    return not violations, violations


def evaluate_record(
    record: StructureRecord,
    ratio_limits: Sequence[Tuple[str, str, float]],
    exclude_counts: Sequence[int],
    composition_targets: Dict[Tuple[Tuple[str, int], ...], str],
) -> Tuple[bool, List[str]]:
    """
    返回 (should_move, reasons)
    """
    if should_skip_by_element_count(record, exclude_counts):
        return False, []

    is_valid, violation_msgs = check_ratio_constraints(record, ratio_limits)
    reasons: List[str] = violation_msgs.copy()

    if composition_targets:
        target = composition_targets.get(record.canonical_key)
        if target:
            reasons.append(f"命中配比 {target}")

    should_move = bool(reasons)

    return should_move, reasons


# -----------------------------------------------------------------------------
# 输出与文件操作
# -----------------------------------------------------------------------------


def print_stats(
    records: Sequence[StructureRecord],
    ratio_limits: Sequence[Tuple[str, str, float]],
    composition_targets: Dict[Tuple[Tuple[str, int], ...], str],
    exclude_counts: Sequence[int],
    output_path: Optional[Path],
) -> None:
    total = len(records)
    composition_counter = Counter(r.signature for r in records)

    print("=" * 70)
    print("配比统计")
    print("=" * 70)
    print(f"结构总数：{total}")
    if ratio_limits:
        print("比例约束：")
        for elem1, elem2, max_ratio in ratio_limits:
            print(f"  - {elem1}:{elem2} <= {max_ratio}:1")
    if composition_targets:
        print("按配比剔除：")
        for display in composition_targets.values():
            print(f"  - {display}")
    if exclude_counts:
        exclude_str = ", ".join(f"{n} 元化合物" for n in sorted(exclude_counts))
        print(f"排除元素数量：{exclude_str}")
    print()

    if ratio_limits or composition_targets:
        exclude_counts = exclude_counts or []
        positives = 0
        negatives = 0
        for record in records:
            should_move, _ = evaluate_record(record, ratio_limits, exclude_counts, composition_targets)
            if should_move:
                positives += 1
            else:
                negatives += 1
        print(f"需剔除：{positives} ({positives * 100 / total:.1f}%)")
        print(f"保留：  {negatives} ({negatives * 100 / total:.1f}%)\n")

    print("配比分布（前 15 项）：")
    for signature, count in composition_counter.most_common(15):
        pct = count * 100 / total if total else 0.0
        print(f"  {signature:<20} {count:>6}  ({pct:>6.2f}%)")
    if len(composition_counter) > 15:
        print(f"  ... 另有 {len(composition_counter) - 15} 种配比")

    if output_path:
        lines = [
            "=" * 70,
            "配比统计报告",
            "=" * 70,
            f"生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}",
            f"结构总数：{total}",
            "",
        ]
        if ratio_limits:
            lines.append("比例约束：")
            for elem1, elem2, max_ratio in ratio_limits:
                lines.append(f"  - {elem1}:{elem2} <= {max_ratio}:1")
            lines.append("")
        if composition_targets:
            lines.append("按配比剔除：")
            for display in composition_targets.values():
                lines.append(f"  - {display}")
            lines.append("")
        lines.append("配比分布：")
        lines.append("-" * 70)
        for signature, count in composition_counter.most_common():
            pct = count * 100 / total if total else 0.0
            lines.append(f"{signature:<20} {count:>6}  ({pct:>6.2f}%)")
        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"\n详情已写入 {output_path}")


def filter_records(
    records: Sequence[StructureRecord],
    ratio_limits: Sequence[Tuple[str, str, float]],
    exclude_counts: Sequence[int],
    composition_targets: Dict[Tuple[Tuple[str, int], ...], str],
) -> None:
    matched: List[StructureRecord] = []
    for record in records:
        should_move, _ = evaluate_record(record, ratio_limits, exclude_counts, composition_targets)
        if not should_move:
            matched.append(record)

    print("=" * 70)
    print(f"满足条件的结构：{len(matched)} / {len(records)}")
    print("=" * 70)
    for rec in matched:
        print(f"{rec.iteration}/{rec.path.name:<40} {rec.signature}")


def confirm(prompt: str, skip: bool) -> bool:
    if skip:
        return True
    try:
        return input(f"{prompt} [y/N]: ").strip().lower() in {"y", "yes"}
    except (KeyboardInterrupt, EOFError):
        print()
        return False


def move_records(
    targets: Sequence[Tuple[StructureRecord, List[str]]],
    trash_dir: Path,
    skip_confirm: bool,
    dry_run: bool,
) -> None:
    if not targets:
        print("没有需要移动的结构。")
        return

    print("将要移动的结构（前 20 个）：")
    for record, reasons in targets[:20]:
        summary = "; ".join(reasons) if reasons else "手动选择"
        print(f"  {record.iteration}/{record.path.name:<35} {record.signature:<15} ({summary})")
    if len(targets) > 20:
        print(f"  ... 另有 {len(targets) - 20} 个结构")
    print()

    if dry_run:
        print("[模拟模式] 不会实际移动文件。")
        return

    task_defs = []
    for record, reasons in targets:
        destination = trash_dir / record.iteration / record.path.name
        task_defs.append((record, reasons, destination))

    if not confirm(f"⚠️  即将移动 {len(task_defs)} 个结构到 {trash_dir}/", skip_confirm):
        print("操作已取消。")
        return

    moved: List[Tuple[Path, Path, List[str]]] = []
    failed: List[Tuple[Path, str]] = []

    workers = available_workers()
    use_parallel = len(task_defs) > 1 and workers > 1

    if use_parallel:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_move_file_task, record.path, destination): (record, reasons, destination)
                for record, reasons, destination in task_defs
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="移动文件", unit="个"):
                record, reasons, destination = futures[future]
                error = future.result()
                if error:
                    failed.append((record.path, error))
                else:
                    moved.append((record.path, destination, reasons))
    else:
        for record, reasons, destination in tqdm(task_defs, desc="移动文件", unit="个"):
            error = _move_file_task(record.path, destination)
            if error:
                failed.append((record.path, error))
            else:
                moved.append((record.path, destination, reasons))

    log_path = trash_dir / "move_log.txt"
    if moved:
        lines = [
            "=" * 70,
            "配比筛选 - 移动日志",
            "=" * 70,
            f"时间：{datetime.now():%Y-%m-%d %H:%M:%S}",
            f"总计：{len(moved)}",
            "",
        ]
        for src, dst, reasons in moved:
            lines.append(f"源文件：{src}")
            lines.append(f"目标：  {dst}")
            if reasons:
                lines.append(f"原因：  {'; '.join(reasons)}")
            lines.append("-" * 70)
        log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"\n成功移动 {len(moved)} 个文件，日志：{log_path}")
    if failed:
        print(f"\n移动失败 {len(failed)} 个文件：")
        for path, msg in failed[:5]:
            print(f"  {path}: {msg}")
        if len(failed) > 5:
            print(f"  ... 另有 {len(failed) - 5} 个失败项")


def restore_trash(trash_dir: Path, skip_confirm: bool) -> None:
    if not trash_dir.exists():
        print(f"垃圾箱目录不存在：{trash_dir}")
        return
    restore_tasks: List[Tuple[Path, Path]] = []
    for it_dir in trash_dir.glob("IT*"):
        if not it_dir.is_dir():
            continue
        for res_path in it_dir.glob("*.res"):
            target = Path(it_dir.name) / res_path.name
            restore_tasks.append((res_path, target))

    if not restore_tasks:
        print("垃圾箱为空，没有需要还原的文件。")
        return

    print(f"将要还原 {len(restore_tasks)} 个文件（显示前 20 个）：")
    for src, dst in restore_tasks[:20]:
        print(f"  {src} → {dst}")
    if len(restore_tasks) > 20:
        print(f"  ... 另有 {len(restore_tasks) - 20} 个文件")
    print()

    if not confirm("⚠️  即将执行还原操作", skip_confirm):
        print("操作已取消。")
        return

    restored = 0
    failed: List[Tuple[Path, str]] = []
    for src, dst in restore_tasks:
        try:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
            restored += 1
        except Exception as exc:  # pragma: no cover - IO 错误
            failed.append((src, str(exc)))

    print(f"已还原 {restored} 个文件。")
    if failed:
        print(f"还原失败 {len(failed)} 个文件：")
        for src, msg in failed[:5]:
            print(f"  {src}: {msg}")
        if len(failed) > 5:
            print(f"  ... 另有 {len(failed) - 5} 个失败项")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于 ASE 的配比筛选工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例：
  python filter_composition_for_PD.py IT0 --stats
  python filter_composition_for_PD.py --ratio-limit Ce:Mg:5 --filter
  python filter_composition_for_PD.py --ratio-limit Ce:Mg:5 -c Ce9Mg8H1 --move-invalid
  python filter_composition_for_PD.py --move-invalid --dry-run
  python filter_composition_for_PD.py --restore --yes
""",
    )

    parser.add_argument(
        "iteration",
        nargs="?",
        help="迭代编号，例如 0 或 IT0；缺省时扫描所有 IT* 目录。",
    )
    parser.add_argument(
        "--ratio-limit",
        action="append",
        metavar="E1:E2:R",
        help="元素比例约束，例如 Ce:Mg:5 表示两者比例不超过 5:1。",
    )
    parser.add_argument(
        "-c",
        "--move-composition",
        nargs="+",
        action="extend",
        metavar="SIGNATURE",
        help="指定需直接剔除的配比签名，可重复使用。",
    )
    parser.add_argument(
        "--exclude-elements",
        type=int,
        nargs="+",
        metavar="N",
        help="指定不参与剔除的元素数量，例如 '1 2' 表示单质和二元化合物跳过。",
    )

    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--stats", action="store_true", help="统计模式")
    mode_group.add_argument("--filter", action="store_true", help="筛选模式（列出满足条件的结构）")
    mode_group.add_argument("--move-invalid", action="store_true", help="移动模式")
    mode_group.add_argument("--restore", action="store_true", help="从垃圾箱还原结构")

    parser.add_argument("--trash-dir", default="trash_bin", help="垃圾箱目录（默认：trash_bin）")
    parser.add_argument("--output", help="统计模式报告输出路径")
    parser.add_argument("--dry-run", action="store_true", help="移动模式下仅展示，不实际移动")
    parser.add_argument("--yes", "-y", action="store_true", help="跳过所有确认提示")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    root = Path.cwd()
    if root.name != "PD":
        print("警告：建议在 PD 目录内运行。当前目录：", root, file=sys.stderr)

    iteration = None
    if args.iteration:
        try:
            iteration = parse_iteration_label(args.iteration)
        except ValueError as exc:
            parser.error(str(exc))

    it_dirs = discover_it_dirs(root, iteration)
    if not it_dirs:
        parser.error("未找到任何 IT* 目录")

    ratio_limits: List[Tuple[str, str, float]] = []
    if args.ratio_limit:
        try:
            ratio_limits = [parse_ratio_limit(item) for item in args.ratio_limit]
        except ValueError as exc:
            parser.error(str(exc))

    composition_targets: Dict[Tuple[Tuple[str, int], ...], str] = {}
    if args.move_composition:
        for signature in args.move_composition:
            key = parse_signature_to_key(signature)
            if key:
                composition_targets[key] = signature.strip()

    exclude_counts = args.exclude_elements or []
    trash_dir = root / args.trash_dir

    if args.restore:
        restore_trash(trash_dir, args.yes)
        return 0

    # 需要结构数据的模式
    records, warnings = load_structures(it_dirs)
    if warnings:
        print("以下文件解析失败：")
        for msg in warnings:
            print(f"  {msg}")
        print()

    if not records:
        print("未成功解析任何结构。")
        return 1

    if args.stats or (not args.filter and not args.move_invalid):
        output_path = Path(args.output) if args.output else None
        print_stats(records, ratio_limits, composition_targets, exclude_counts, output_path)
        if not (args.filter or args.move_invalid):
            return 0

    if args.filter:
        filter_records(records, ratio_limits, exclude_counts, composition_targets)
        if not args.move_invalid:
            return 0

    if args.move_invalid:
        targets: List[Tuple[StructureRecord, List[str]]] = []
        for record in records:
            should_move, reasons = evaluate_record(
                record, ratio_limits, exclude_counts, composition_targets
            )
            if should_move:
                targets.append((record, reasons))

        move_records(targets, trash_dir, args.yes, args.dry_run)

    return 0


if __name__ == "__main__":
    sys.exit(main())
