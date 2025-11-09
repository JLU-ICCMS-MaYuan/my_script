#!/usr/bin/env python3
"""
Utilities for summarising LAMMPS relaxation jobs inside the RELAX directory.

请直接在 RELAX 目录中运行本脚本，例如：
    python analyze_relax.py            # 遍历当前目录下所有 IT* 子目录
    python analyze_relax.py 0 3 5      # 仅检查 IT0、IT3、IT5
    python analyze_relax.py IT2 IT4    # IT 前缀可选
"""

from __future__ import annotations

import argparse
import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback when tqdm not installed
    def tqdm(iterable, *args, **kwargs):  # type: ignore
        return iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="统计当前 RELAX 目录内各 IT* 代的结构优化情况。"
    )
    parser.add_argument(
        "iterations",
        nargs="*",
        help="指定要检查的代号（可以写成 3 或 IT3）。缺省时遍历所有 IT* 目录。",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="并行分析线程数；缺省或 <=0 时自动按 CPU 核心数设置。",
    )
    return parser.parse_args()


def normalise_iteration(token: str) -> int:
    match = re.fullmatch(r"(?:IT)?(\d+)", token.strip(), re.IGNORECASE)
    if not match:
        raise ValueError(f"无法解析的 IT 编号: {token!r}")
    return int(match.group(1))


def discover_iterations(relax_root: Path, tokens: Iterable[str]) -> List[Tuple[int, Path]]:
    if tokens:
        wanted = []
        for tok in tokens:
            idx = normalise_iteration(tok)
            it_dir = relax_root / f"IT{idx}"
            wanted.append((idx, it_dir))
    else:
        wanted = []
        for child in relax_root.iterdir():
            if child.is_dir():
                match = re.fullmatch(r"IT(\d+)", child.name)
                if match:
                    wanted.append((int(match.group(1)), child))
    wanted.sort(key=lambda item: item[0])
    return wanted


def load_lammps_targets(it_dir: Path) -> OrderedDict[str, List[str]]:
    targets: OrderedDict[str, List[str]] = OrderedDict()
    group_files = sorted(it_dir.glob("group_*"))
    for gf in group_files:
        with gf.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                name = Path(line).stem
                targets.setdefault(name, []).append(gf.name)
    return targets


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def classify_structure(it_dir: Path, struct_name: str) -> Tuple[str, str | None]:
    struct_dir = it_dir / struct_name
    if not struct_dir.exists():
        return "not_started", "missing directory"

    log_path = struct_dir / "log.lammps"
    if not log_path.exists():
        return "not_started", "missing log.lammps"

    log_text = read_text(log_path)
    has_opt_end = "<OPTEND>" in log_text
    has_opt_fail = "<OPTFAIL>" in log_text

    out_res = struct_dir / f"{struct_name}-out.res"
    limit_res = struct_dir / f"{struct_name}-limitbonds.res"

    if out_res.exists() and has_opt_end and not has_opt_fail:
        return "success", None

    if has_opt_fail:
        reason = "log contains <OPTFAIL>"
    elif limit_res.exists():
        reason = "limitbonds file"
    elif not has_opt_end:
        reason = "log missing <OPTEND>"
    elif not out_res.exists():
        reason = "missing -out.res"
    else:
        reason = "unknown"
    return "failed", reason


def detect_iteration_mode(it_dir: Path) -> str:
    has_bfgs = any(it_dir.glob("relax_group_*"))
    if has_bfgs:
        return "bfgs"
    return "lammps"


def parse_bfgs_group_file(path: Path) -> Tuple[List[str], Dict[str, Tuple[str, str | None]]]:
    label_pattern = re.compile(r"\s*Label\s+(\d+)\s+(\S+)")
    frame_pattern = re.compile(r"#optimization frame\s+(\d+)")

    labels_map: Dict[int, str] = {}
    statuses: Dict[str, Tuple[str, str | None]] = {}
    seen_frames: Set[int] = set()

    current_frame: int | None = None
    frame_success = False
    frame_failure_reason: str | None = None

    def finalize_frame() -> None:
        nonlocal current_frame, frame_success, frame_failure_reason
        if current_frame is None:
            return
        name = labels_map.get(current_frame)
        if name is None:
            current_frame = None
            frame_success = False
            frame_failure_reason = None
            return
        seen_frames.add(current_frame)
        if frame_success:
            statuses[name] = ("success", None)
        else:
            reason = frame_failure_reason or "BFGS 日志缺少成功标记"
            statuses[name] = ("failed", reason)
        current_frame = None
        frame_success = False
        frame_failure_reason = None

    failure_keywords = ("error", "fail", "diverg", "nan", "abort")

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            label_match = label_pattern.match(line)
            if label_match:
                idx = int(label_match.group(1)) - 1
                name = label_match.group(2)
                labels_map[idx] = name
                continue

            frame_match = frame_pattern.match(line)
            if frame_match:
                finalize_frame()
                current_frame = int(frame_match.group(1))
                frame_success = False
                frame_failure_reason = None
                continue

            stripped = line.strip()
            if not stripped:
                continue

            lowered = stripped.lower()
            if "relaxation is converged!" in lowered:
                frame_success = True
                frame_failure_reason = None
                continue

            if current_frame is not None:
                if any(keyword in lowered for keyword in failure_keywords):
                    frame_failure_reason = stripped
                elif "relaxation is not converged yet!" in lowered:
                    frame_failure_reason = "Relaxation 未达到收敛判据"

    finalize_frame()

    if labels_map:
        max_idx = max(labels_map)
        ordered_labels: List[str] = ["" for _ in range(max_idx + 1)]
        for idx, name in labels_map.items():
            ordered_labels[idx] = name
    else:
        ordered_labels = []

    for frame_idx, name in enumerate(ordered_labels):
        if not name:
            continue
        if name not in statuses:
            if frame_idx in seen_frames:
                statuses[name] = ("failed", "BFGS 日志缺少成功标记")
            else:
                statuses[name] = ("not_started", "BFGS 日志缺少该结构的优化记录")

    return ordered_labels, statuses


def collect_bfgs_data(
    it_dir: Path, requested_workers: int | None
) -> Tuple[OrderedDict[str, List[str]], Dict[str, Tuple[str, str | None]]]:
    group_files = sorted(it_dir.glob("relax_group_*"))
    if not group_files:
        return OrderedDict(), {}

    def parse_file(path: Path) -> Tuple[str, List[str], Dict[str, Tuple[str, str | None]]]:
        labels, statuses = parse_bfgs_group_file(path)
        return path.name, labels, statuses

    workers = resolve_worker_count(len(group_files), requested_workers)
    results: List[Tuple[str, List[str], Dict[str, Tuple[str, str | None]]]] = []
    if workers <= 1:
        for gf in group_files:
            results.append(parse_file(gf))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for item in executor.map(parse_file, group_files):
                results.append(item)

    targets: OrderedDict[str, List[str]] = OrderedDict()
    combined_statuses: Dict[str, Tuple[str, str | None]] = {}
    for file_name, labels, status_map in results:
        for name in labels:
            if not name:
                continue
            targets.setdefault(name, []).append(file_name)
        combined_statuses.update(status_map)

    return targets, combined_statuses


def classify_lammps_targets(
    it_dir: Path, target_names: List[str], workers: int
) -> Dict[str, Tuple[str, str | None]]:
    def classify_with_name(name: str) -> Tuple[str, Tuple[str, str | None]]:
        status, reason = classify_structure(it_dir, name)
        return name, (status, reason)

    results: Dict[str, Tuple[str, str | None]] = {}
    if workers <= 1:
        for struct_name in target_names:
            name, outcome = classify_with_name(struct_name)
            results[name] = outcome
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for name, outcome in executor.map(classify_with_name, target_names):
                results[name] = outcome
    return results


def resolve_worker_count(total: int, requested: int | None) -> int:
    if total <= 0:
        raise ValueError("total must be positive for worker allocation")
    if requested is not None and requested > 0:
        return min(requested, total)
    cpu_count = os.cpu_count() or 1
    return max(1, min(cpu_count, total))


def write_iteration_log(
    it_dir: Path,
    it_index: int,
    stats: Dict[str, List[Tuple[str, str | None]]],
    duplicates: Dict[str, List[str]],
    mode: str,
) -> Path:
    log_path = it_dir / f"IT{it_index}_relax_summary.log"
    lines: List[str] = []
    total = sum(len(items) for items in stats.values())
    lines.append(f"Iteration IT{it_index}")
    lines.append(f"Optimization mode: {mode.upper()}")
    lines.append(f"Total structures: {total}")
    lines.append(f"Success: {len(stats['success'])}")
    lines.append(f"Failed: {len(stats['failed'])}")
    lines.append(f"Not started: {len(stats['not_started'])}")
    if duplicates:
        lines.append("")
        lines.append("Duplicates detected (structure listed in multiple source files):")
        for name in sorted(duplicates):
            sources = ", ".join(sorted(duplicates[name]))
            lines.append(f"  {name} -> {sources}")
    category_titles = {
        "success": "[Success]",
        "failed": "[Failed]",
        "not_started": "[Not started]",
    }
    for category in ("success", "failed", "not_started"):
        items = stats[category]
        lines.append("")
        lines.append(category_titles[category])
        if not items:
            lines.append("  (none)")
            continue
        for name, reason in items:
            if reason:
                lines.append(f"  {name} -- {reason}")
            else:
                lines.append(f"  {name}")
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return log_path


def write_category_files(
    it_dir: Path,
    it_index: int,
    stats: Dict[str, List[Tuple[str, str | None]]],
) -> Tuple[Path, Path]:
    success_path = it_dir / f"IT{it_index}_relax_success.txt"
    failed_path = it_dir / f"IT{it_index}_relax_failed.txt"

    success_entries = [name for name, _ in stats["success"]] or ["(none)"]
    failed_entries = [
        f"{name} -- {reason}" if reason else name for name, reason in stats["failed"]
    ] or ["(none)"]

    success_path.write_text("\n".join(success_entries) + "\n", encoding="utf-8")
    failed_path.write_text("\n".join(failed_entries) + "\n", encoding="utf-8")
    return success_path, failed_path


def main() -> None:
    args = parse_args()
    relax_root = Path.cwd()
    if relax_root.name != "RELAX":
        raise SystemExit(
            f"当前目录应为 'RELAX'，但检测到: {relax_root}. 请切换到 RELAX 目录后重试。"
        )

    iterations = discover_iterations(relax_root, args.iterations)
    if not iterations:
        raise SystemExit("未在指定目录下找到任何 IT* 子目录。")

    for idx, it_dir in iterations:
        mode = detect_iteration_mode(it_dir)

        if mode == "bfgs":
            targets, status_map = collect_bfgs_data(it_dir, args.workers)
            if not targets:
                print(f"IT{idx}: 未找到任何 relax_group_* 文件，跳过。")
                continue
        else:
            targets = load_lammps_targets(it_dir)
            if not targets:
                print(f"IT{idx}: 未找到任何 group_* 文件，跳过。")
                continue
            target_names = list(targets)
            workers = resolve_worker_count(len(target_names), args.workers)
            status_map = classify_lammps_targets(it_dir, target_names, workers)

        duplicates = {name: sources for name, sources in targets.items() if len(sources) > 1}

        stats: Dict[str, List[Tuple[str, str | None]]] = {
            "success": [],
            "failed": [],
            "not_started": [],
        }

        target_names = list(targets)
        for name, (status, reason) in tqdm(
            ((struct_name, status_map.get(struct_name, ("not_started", "未知状态")))
             for struct_name in target_names),
            total=len(target_names),
            desc=f"IT{idx}",
            unit="结构",
        ):
            stats[status].append((name, reason))

        for category in stats:
            stats[category].sort(key=lambda item: item[0])

        log_file = write_iteration_log(it_dir, idx, stats, duplicates, mode)
        success_file, failed_file = write_category_files(it_dir, idx, stats)

        total = sum(len(items) for items in stats.values())
        print(
            f"IT{idx}({mode.upper()}): total={total} success={len(stats['success'])} "
            f"failed={len(stats['failed'])} not_started={len(stats['not_started'])} "
            f"-> {log_file}, success->{success_file}, failed->{failed_file}"
        )


if __name__ == "__main__":
    main()
