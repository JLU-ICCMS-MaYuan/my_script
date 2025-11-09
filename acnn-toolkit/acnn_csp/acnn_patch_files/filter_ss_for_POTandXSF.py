#!/usr/bin/env python3
"""过滤 POT/IT*/ss 中超过阈值的结构，并列出 XSF 对应路径。"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

# 浮点匹配正则
FLOAT_RE = r"[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?"
LINE_RE = re.compile(
    rf"^(?P<name>\S+\.xsf)\s+\d+\s+(?P<e>{FLOAT_RE})\s+(?P<f>{FLOAT_RE})\s+(?P<v>{FLOAT_RE})"
)


@dataclass(frozen=True)
class Entry:
    name: str
    e: float
    f: float
    v: float

    def exceeds(self, e_lim: Optional[float], f_lim: Optional[float], v_lim: Optional[float]) -> bool:
        if e_lim is not None and self.e > e_lim:
            return True
        if f_lim is not None and self.f > f_lim:
            return True
        if v_lim is not None and self.v > v_lim:
            return True
        return False

    def exceeded_fields(
        self, e_lim: Optional[float], f_lim: Optional[float], v_lim: Optional[float]
    ) -> Sequence[str]:
        fields: List[str] = []
        if e_lim is not None and self.e > e_lim:
            fields.append("e")
        if f_lim is not None and self.f > f_lim:
            fields.append("f")
        if v_lim is not None and self.v > v_lim:
            fields.append("v")
        return fields


def parse_iteration_label(label: str) -> str:
    value = label.strip()
    if not value:
        raise ValueError("迭代编号不能为空")
    value_upper = value.upper()
    if value_upper.startswith("IT"):
        return value_upper
    return f"IT{value}"


def load_entries(ss_path: Path) -> List[Entry]:
    entries: List[Entry] = []
    with ss_path.open("r", encoding="ascii", errors="ignore") as fh:
        for line in fh:
            match = LINE_RE.match(line)
            if not match:
                continue
            try:
                entry = Entry(
                    name=match.group("name"),
                    e=float(match.group("e")),
                    f=float(match.group("f")),
                    v=float(match.group("v")),
                )
            except ValueError:
                continue
            entries.append(entry)
    return entries


def index_xsf_paths(xsf_root: Path) -> dict[str, List[Path]]:
    file_map: dict[str, List[Path]] = {}
    for it_dir in xsf_root.glob("IT*"):
        if not it_dir.is_dir():
            continue
        for xsf_file in it_dir.rglob("*.xsf"):
            if not xsf_file.is_file():
                continue
            name = xsf_file.name
            resolved_file = xsf_file.resolve()
            file_map.setdefault(name, []).append(resolved_file)
    return file_map


def existing_paths_for(entry: Entry, file_map: dict[str, List[Path]], xsf_root: Path) -> List[Path]:
    paths: List[Path] = []
    seen = set()
    for candidate in file_map.get(entry.name, []):
        if candidate not in seen:
            paths.append(candidate)
            seen.add(candidate)

    if not paths:
        # 兜底：尝试在所有 IT 目录直接定位
        for xsf_file in xsf_root.glob(f"IT*/{entry.name}"):
            resolved = xsf_file.resolve()
            if resolved not in seen and resolved.is_file():
                paths.append(resolved)
                seen.add(resolved)
    return paths


def relpath(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return os.path.relpath(path, start=base)


def delete_xsf_files(paths: Sequence[Path]) -> Tuple[int, int]:
    files = list(paths)
    total = len(files)
    if total == 0:
        return 0, 0
    removed = 0
    missing = 0
    bar_width = 40
    for idx, xsf_file in enumerate(files, start=1):
        try:
            xsf_file.unlink()
        except FileNotFoundError:
            missing += 1
        except IsADirectoryError:
            missing += 1
        else:
            removed += 1
        filled = int(bar_width * idx / total)
        bar = "#" * filled + "." * (bar_width - filled)
        sys.stdout.write(f"\r删除进度：[{bar}] {idx}/{total}")
        sys.stdout.flush()
    sys.stdout.write("\n")
    return removed, missing


def move_dt_xsf(it_dir: Path, entries: Iterable[Entry]) -> Tuple[int, int, int]:
    dt_dir = it_dir / "DT"
    if not dt_dir.is_dir():
        return 0, 0, 0
    trashbin = it_dir / "trashbin"
    trashbin.mkdir(exist_ok=True)
    moved = 0
    skipped = 0
    missing = 0
    handled: set[str] = set()
    for entry in entries:
        if entry.name in handled:
            continue
        handled.add(entry.name)
        source = dt_dir / entry.name
        if not source.exists():
            missing += 1
            continue
        dest = trashbin / source.name
        if dest.exists():
            skipped += 1
            continue
        source.rename(dest)
        moved += 1
    return moved, skipped, missing


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="在 POT 目录中运行，指定 IT 编号及阈值，筛选超限结构。"
    )
    parser.add_argument(
        "iteration",
        help="迭代编号，例如 0 或 IT0。",
    )
    parser.add_argument(
        "--xsf-root",
        default="../XSF",
        help="XSF 根目录（相对或绝对路径，默认 ../XSF）。",
    )
    parser.add_argument(
        "--e-limit",
        type=float,
        help="e/atom 上限。",
    )
    parser.add_argument(
        "--f-limit",
        type=float,
        help="f/atom 上限。",
    )
    parser.add_argument(
        "--v-limit",
        type=float,
        help="v/atom 上限。",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="结果输出文件（默认 IT*/exceeded_structures.txt）。",
    )
    parser.add_argument(
        "-d",
        "--delete",
        action="store_true",
        help="删除匹配的 .xsf 文件（删除过程中显示进度）。",
    )
    parser.add_argument(
        "--move-dt",
        action="store_true",
        help="将超限结构在 DT 目录中的 .xsf 文件移动到同迭代的 trashbin 目录。",
    )
    args = parser.parse_args(argv)
    if args.e_limit is None and args.f_limit is None and args.v_limit is None:
        parser.error("至少需要提供 --e-limit、--f-limit、--v-limit 中的一个。")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    pot_dir = Path.cwd()
    if pot_dir.name != "POT":
        print("警告：当前路径不是 POT 目录，建议切换到 POT 后再运行。", file=sys.stderr)

    try:
        iteration = parse_iteration_label(args.iteration)
    except ValueError as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 1

    it_dir = pot_dir / iteration
    if not it_dir.is_dir():
        print(f"错误：未找到迭代目录 {it_dir}", file=sys.stderr)
        return 1

    ss_path = it_dir / "ss"
    if not ss_path.is_file():
        print(f"错误：未找到 ss 文件 {ss_path}", file=sys.stderr)
        return 1

    entries = load_entries(ss_path)
    if not entries:
        print(f"警告：未能在 {ss_path} 中解析到有效数据。", file=sys.stderr)

    xsf_root = Path(args.xsf_root)
    if not xsf_root.is_absolute():
        xsf_root = (pot_dir / xsf_root).resolve()
    if not xsf_root.is_dir():
        print(f"错误：未找到 XSF 根目录 {xsf_root}", file=sys.stderr)
        return 1

    file_map = index_xsf_paths(xsf_root)

    matched: List[Tuple[Entry, List[Path]]] = []
    for entry in entries:
        if not entry.exceeds(args.e_limit, args.f_limit, args.v_limit):
            continue
        paths = existing_paths_for(entry, file_map, xsf_root)
        matched.append((entry, paths))

    if not matched:
        print("没有结构超过设定阈值。")
        return 0

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = it_dir / "exceeded_structures.txt"
    if not output_path.is_absolute():
        output_path = (pot_dir / output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    for entry, paths in matched:
        reasons = ",".join(entry.exceeded_fields(args.e_limit, args.f_limit, args.v_limit))
        if paths:
            for path in paths:
                path_str = relpath(path, pot_dir)
                lines.append(
                    f"{path_str:<110}  {reasons:<8}  {entry.e:.6g}  {entry.f:.6g}  {entry.v:.6g}"
                )
        else:
            path_str = f"# 未找到对应路径 {entry.name}"
            lines.append(
                f"{path_str:<110}  {reasons:<8}  {entry.e:.6g}  {entry.f:.6g}  {entry.v:.6g}"
            )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"已写入 {len(lines)} 条记录：{relpath(output_path, pot_dir)}")

    if args.delete:
        target_files = sorted(
            {
                path
                for _, paths in matched
                for path in paths
                if path.suffix.lower() == ".xsf"
            }
        )
        removed_count, missing_count = delete_xsf_files(target_files)
        if removed_count:
            if missing_count:
                print(f"已删除 {removed_count} 个 XSF 文件，另有 {missing_count} 个未找到。")
            else:
                print(f"已删除 {removed_count} 个 XSF 文件。")
        else:
            if missing_count:
                print(f"未找到 {missing_count} 个待删除的 XSF 文件。")
            else:
                print("未找到可删除的 XSF 文件。")

    missing = sum(1 for _, paths in matched if not paths)
    if missing:
        print(f"注意：{missing} 个结构在 XSF 中未找到对应路径。")

    if args.move_dt:
        matched_entries = [entry for entry, _ in matched]
        moved, skipped, missing = move_dt_xsf(it_dir, matched_entries)
        if moved or skipped or missing:
            print(
                "DT 中超限结构已处理："
                f"移动 {moved} 个，跳过 {skipped} 个（trashbin 已存在同名文件），"
                f"未找到 {missing} 个。"
            )
        else:
            print("DT 目录不存在或无可移动的超限结构 .xsf 文件。")

    return 0


if __name__ == "__main__":
    sys.exit(main())
