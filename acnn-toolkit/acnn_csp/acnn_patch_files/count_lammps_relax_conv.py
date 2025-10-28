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
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="统计当前 RELAX 目录内各 IT* 代的结构优化情况。"
    )
    parser.add_argument(
        "iterations",
        nargs="*",
        help="指定要检查的代号（可以写成 3 或 IT3）。缺省时遍历所有 IT* 目录。",
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


def load_targets(it_dir: Path) -> OrderedDict[str, List[str]]:
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


def write_iteration_log(
    it_dir: Path,
    it_index: int,
    stats: Dict[str, List[Tuple[str, str | None]]],
    duplicates: Dict[str, List[str]],
) -> Path:
    log_path = it_dir / f"IT{it_index}_relax_summary.log"
    lines: List[str] = []
    total = sum(len(items) for items in stats.values())
    lines.append(f"Iteration IT{it_index}")
    lines.append(f"Total structures: {total}")
    lines.append(f"Success: {len(stats['success'])}")
    lines.append(f"Failed: {len(stats['failed'])}")
    lines.append(f"Not started: {len(stats['not_started'])}")
    if duplicates:
        lines.append("")
        lines.append("Duplicates detected (structure listed in multiple group_* files):")
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
        targets = load_targets(it_dir)
        if not targets:
            print(f"IT{idx}: 未找到任何 group_* 文件，跳过。")
            continue

        duplicates = {name: sources for name, sources in targets.items() if len(sources) > 1}

        stats: Dict[str, List[Tuple[str, str | None]]] = {
            "success": [],
            "failed": [],
            "not_started": [],
        }

        for name in targets:
            status, reason = classify_structure(it_dir, name)
            stats[status].append((name, reason))

        for category in stats:
            stats[category].sort(key=lambda item: item[0])

        log_file = write_iteration_log(it_dir, idx, stats, duplicates)

        total = sum(len(items) for items in stats.values())
        print(
            f"IT{idx}: total={total} success={len(stats['success'])} "
            f"failed={len(stats['failed'])} not_started={len(stats['not_started'])} "
            f"-> {log_file}"
        )


if __name__ == "__main__":
    main()
