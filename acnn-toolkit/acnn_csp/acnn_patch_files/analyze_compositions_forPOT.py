#!/usr/bin/env python3
"""Analyze ternary compositions and generate interactive reports."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import plotly.express as px
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass
class StructureInfo:
    name: str
    counts: Tuple[int, int, int]


@dataclass
class ErrorRecord:
    name: str
    energy: float
    force: float
    virial: float
    composition: str = ""
    over_energy: bool = False
    over_force: bool = False
    over_virial: bool = False


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="统计配比并生成三角图")
    parser.add_argument("iteration", help="迭代编号，例如 IT0")
    parser.add_argument("--ss-path", help="自定义 ss 文件路径")
    parser.add_argument("--output-dir", default=None, help="输出目录（默认使用迭代目录）")
    parser.add_argument(
        "--html",
        default="composition_summary.html",
        help="主配比图 HTML 文件名（默认 composition_summary.html）",
    )
    parser.add_argument(
        "--high-html",
        default="beyond_rmse.html",
        help="超 RMSE 图 HTML 文件名（默认 beyond_rmse.html）",
    )
    return parser.parse_args()


def log_step(message: str) -> None:
    print(f"[进度] {message}")


# ---------------------------------------------------------------------------
# Structure parsing helpers
# ---------------------------------------------------------------------------


def read_xsf_composition(path: Path) -> Counter[str]:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()
    for idx, line in enumerate(lines):
        if line.strip().upper().startswith("PRIMCOORD") and idx + 1 < len(lines):
            n_atoms = int(lines[idx + 1].split()[0])
            composition: Counter[str] = Counter()
            for atom_line in lines[idx + 2 : idx + 2 + n_atoms]:
                parts = atom_line.split()
                if parts:
                    composition[parts[0]] += 1
            if sum(composition.values()) != n_atoms:
                raise ValueError("PRIMCOORD 段长度不足")
            return composition
    raise ValueError(f"未找到 PRIMCOORD 段：{path}")


def reduce_tuple(counts: Tuple[int, ...]) -> Tuple[int, ...]:
    gcd_val = 0
    for value in counts:
        gcd_val = math.gcd(gcd_val, value)
    gcd_val = gcd_val or 1
    return tuple(value // gcd_val for value in counts)


def tuple_to_formula(elements: Sequence[str], counts: Tuple[int, ...]) -> str:
    return "".join(f"{elem}{count}" for elem, count in zip(elements, counts))


def load_structures(dt_dir: Path) -> Tuple[Dict[str, StructureInfo], List[str]]:
    compositions: Dict[str, Counter[str]] = {}
    elements_set: set[str] = set()
    for xsf in sorted(dt_dir.glob("*.xsf")):
        try:
            comp = read_xsf_composition(xsf)
        except ValueError as exc:
            print(f"跳过解析失败的文件：{xsf.name} ({exc})")
            continue
        compositions[xsf.name] = comp
        elements_set.update(comp.keys())

    if not compositions:
        raise RuntimeError("DT 目录中没有可解析的 .xsf 文件")

    elements = sorted(elements_set)
    if len(elements) != 3:
        raise RuntimeError(f"检测到 {len(elements)} 种元素，脚本仅支持三元体系")

    structures: Dict[str, StructureInfo] = {}
    for name, comp in compositions.items():
        counts = tuple(comp.get(elem, 0) for elem in elements)
        structures[name] = StructureInfo(name=name, counts=counts)

    return structures, elements


# ---------------------------------------------------------------------------
# ss parsing
# ---------------------------------------------------------------------------


def parse_ss_file(path: Path) -> Tuple[List[ErrorRecord], Dict[str, float]]:
    rmse_pattern = re.compile(r"rmse([EFV])\s*:?\s*([0-9.eE+-]+)")
    entry_pattern = re.compile(
        r"([\w.-]+\.xsf)\s+\d+\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)"
    )

    thresholds: Dict[str, float] = {}
    records: List[ErrorRecord] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            for match in rmse_pattern.finditer(line):
                mapping = {"E": "energy", "F": "force", "V": "virial"}
                thresholds.setdefault(mapping[match.group(1)], float(match.group(2)))

            match = entry_pattern.search(line)
            if not match:
                continue
            name, err_e, err_f, err_v = match.groups()
            records.append(
                ErrorRecord(
                    name=name,
                    energy=float(err_e),
                    force=float(err_f),
                    virial=float(err_v),
                )
            )

    if not thresholds:
        raise RuntimeError("ss 文件中未找到 rmse 信息")

    return records, thresholds


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def write_sorted_errors(
    records: List[ErrorRecord],
    attr: str,
    path: Path,
) -> None:
    sorted_records = sorted(records, key=lambda item: getattr(item, attr), reverse=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("结构\t配比\t误差\n")
        for record in sorted_records:
            handle.write(
                f"{record.name:<66}\t{record.composition}\t{getattr(record, attr):.6e}\n"
            )


def compute_rmse(values: Iterable[float]) -> float:
    seq = list(values)
    if not seq:
        return float("nan")
    return math.sqrt(sum(v * v for v in seq) / len(seq))


def summarize_errors(
    structures: Dict[str, StructureInfo],
    entries: List[ErrorRecord],
    elements: Sequence[str],
    thresholds: Dict[str, float],
):
    ratio_map = {info.name: info.counts for info in structures.values()}
    composition_counter: Dict[Tuple[int, ...], int] = defaultdict(int)
    ratio_errors: Dict[Tuple[int, ...], Dict[str, List[float]]] = defaultdict(lambda: defaultdict(list))
    high_ratio_counts: Dict[Tuple[int, ...], int] = defaultdict(int)
    metric_counts = {
        "energy": defaultdict(int),
        "force": defaultdict(int),
        "virial": defaultdict(int),
    }
    high_structures: List[ErrorRecord] = []

    for info in structures.values():
        composition_counter[info.counts] += 1

    for entry in entries:
        counts = ratio_map.get(entry.name)
        if counts is None:
            continue
        entry.composition = tuple_to_formula(elements, counts)
        ratio_errors[counts]["energy"].append(entry.energy)
        ratio_errors[counts]["force"].append(entry.force)
        ratio_errors[counts]["virial"].append(entry.virial)

        if thresholds.get("energy") and entry.energy > thresholds["energy"]:
            entry.over_energy = True
            metric_counts["energy"][counts] += 1
        if thresholds.get("force") and entry.force > thresholds["force"]:
            entry.over_force = True
            metric_counts["force"][counts] += 1
        if thresholds.get("virial") and entry.virial > thresholds["virial"]:
            entry.over_virial = True
            metric_counts["virial"][counts] += 1

        if entry.over_energy or entry.over_force or entry.over_virial:
            high_ratio_counts[counts] += 1
            high_structures.append(entry)

    return composition_counter, ratio_errors, high_ratio_counts, metric_counts, high_structures


def write_composition_stats(
    path: Path,
    elements: Sequence[str],
    composition_counter: Dict[Tuple[int, ...], int],
    ratio_errors: Dict[Tuple[int, ...], Dict[str, List[float]]],
) -> None:
    rows: List[Tuple[str, int, float, float, float]] = []
    for ratio, count in composition_counter.items():
        err = ratio_errors.get(ratio, {})
        rmse_e = compute_rmse(err.get("energy", []))
        rmse_f = compute_rmse(err.get("force", []))
        rmse_v = compute_rmse(err.get("virial", []))
        rows.append((tuple_to_formula(elements, ratio), count, rmse_e, rmse_f, rmse_v))
    rows.sort(key=lambda item: item[1], reverse=True)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("配比\t结构数\tRMSE_E\tRMSE_F\tRMSE_V\n")
        for row in rows:
            handle.write(
                f"{row[0]}\t{row[1]}\t{row[2]:.6e}\t{row[3]:.6e}\t{row[4]:.6e}\n"
            )


def write_high_error_structures(records: List[ErrorRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("结构\t配比\tenergy\tforce\tvirial\n")
        for record in records:
            handle.write(
                f"{record.name:<66}\t{record.composition}\t{record.energy:.6e}\t{record.force:.6e}\t{record.virial:.6e}\n"
            )


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def ratio_label(value: float, elem1: str, elem2: str) -> str:
    frac = Fraction(value).limit_denominator(5)
    numerator = frac.numerator
    denominator = frac.denominator
    label = []
    if numerator:
        label.append(f"{elem1}{numerator}")
    if denominator - numerator:
        label.append(f"{elem2}{denominator - numerator}")
    return "".join(label)


def annotate_edges(fig: go.Figure, elements: Sequence[str]) -> None:
    ticks = [0.2, 0.4, 0.6, 0.8]
    offset = 0.02
    scale = 0.95

    def add_text(a_vals, b_vals, c_vals, texts):
        fig.add_trace(
            go.Scatterternary(
                a=a_vals,
                b=b_vals,
                c=c_vals,
                mode="text",
                text=texts,
                hoverinfo="skip",
                textposition="middle center",
                showlegend=False,
            )
        )

    combos = [
        (elements[1], elements[2], lambda t: (offset, t * scale + offset, (1 - t) * scale + offset)),
        (elements[0], elements[2], lambda t: (t * scale + offset, offset, (1 - t) * scale + offset)),
        (elements[0], elements[1], lambda t: (t * scale + offset, (1 - t) * scale + offset, offset)),
    ]
    for elem_a, elem_b, coord_fn in combos:
        texts, a_vals, b_vals, c_vals = [], [], [], []
        for t in ticks:
            label = ratio_label(t, elem_a, elem_b)
            if not label:
                continue
            a, b, c = coord_fn(t)
            a_vals.append(a)
            b_vals.append(b)
            c_vals.append(c)
            texts.append(label)
        if texts:
            add_text(a_vals, b_vals, c_vals, texts)


def build_dataset(
    ratio_counts: Dict[Tuple[int, ...], int],
    elements: Sequence[str],
) -> List[Dict[str, float]]:
    dataset = []
    for ratio, count in ratio_counts.items():
        total = sum(ratio)
        if not total:
            continue
        dataset.append(
            {
                "a": ratio[0] / total,
                "b": ratio[1] / total,
                "c": ratio[2] / total,
                "count": count,
                "composition": tuple_to_formula(elements, ratio),
            }
        )
    return dataset


def build_ternary_figure(
    dataset: List[Dict[str, float]],
    elements: Sequence[str],
    title: str,
    color_key: str = "count",
    color_title: str = "结构数",
) -> go.Figure | None:
    if not dataset:
        return None
    fig = px.scatter_ternary(
        dataset,
        a="a",
        b="b",
        c="c",
        color=color_key,
        hover_name="composition",
        hover_data={color_key: True},
        color_continuous_scale="Turbo",
        title=title,
    )
    fig.update_traces(marker=dict(size=8))
    fig.update_layout(
        ternary=dict(
            aaxis=dict(title=elements[0]),
            baxis=dict(title=elements[1]),
            caxis=dict(title=elements[2]),
        ),
        coloraxis_colorbar=dict(title=color_title),
    )
    annotate_edges(fig, elements)
    return fig


def export_html(figures: List[go.Figure], path: Path) -> None:
    if not figures:
        return
    html = []
    for idx, fig in enumerate(figures):
        include_js = "cdn" if idx == 0 else False
        html.append(
            fig.to_html(include_plotlyjs=include_js, full_html=False, auto_play=False)
        )
    path.write_text("\n".join(html), encoding="utf-8")


def build_metric_figures(
    elements: Sequence[str],
    metric_name: str,
    count_map: Dict[Tuple[int, ...], int],
    structures: Dict[str, StructureInfo],
    high_structures: List[ErrorRecord],
    flag_attr: str,
) -> List[go.Figure]:
    figures: List[go.Figure] = []
    struct_dataset = []
    for record in high_structures:
        if not getattr(record, flag_attr):
            continue
        info = structures.get(record.name)
        if info is None:
            continue
        total = sum(info.counts)
        if not total:
            continue
        ratio = info.counts
        struct_dataset.append(
            {
                "a": info.counts[0] / total,
                "b": info.counts[1] / total,
                "c": info.counts[2] / total,
                "composition": record.composition,
                "structure": record.name,
                "count": count_map.get(ratio, 0),
            }
        )

    if struct_dataset:
        fig_struct = build_ternary_figure(
            struct_dataset,
            elements,
            f"{metric_name} - 结构",
            color_key="count",
            color_title="结构数",
        )
        if fig_struct:
            figures.append(fig_struct)
    return figures


def parse_iteration_index(label: str) -> int:
    match = re.match(r"(?i)IT(\d+)", label)
    if not match:
        raise ValueError(f"无法解析迭代编号：{label}")
    return int(match.group(1))


def collect_previous_ratios(
    root: Path,
    elements: Sequence[str],
    target_label: str,
) -> set[Tuple[int, ...]]:
    target_idx = parse_iteration_index(target_label)
    ratios: set[Tuple[int, ...]] = set()
    for dir_path in sorted(root.glob("IT*")):
        try:
            idx = parse_iteration_index(dir_path.name)
        except ValueError:
            continue
        if idx >= target_idx:
            continue
        dt_dir = dir_path / "DT"
        if not dt_dir.is_dir():
            continue
        for xsf in dt_dir.glob("*.xsf"):
            try:
                comp = read_xsf_composition(xsf)
            except Exception:
                continue
            counts = tuple(comp.get(elem, 0) for elem in elements)
            ratios.add(counts)
    return ratios


def write_new_compositions(
    path: Path,
    new_infos: List[StructureInfo],
    elements: Sequence[str],
    error_map: Dict[str, ErrorRecord],
) -> None:
    formatter = lambda info: tuple_to_formula(elements, info.counts)
    new_infos = sorted(new_infos, key=lambda info: (formatter(info), info.name))
    with path.open("w", encoding="utf-8") as handle:
        handle.write("结构\t配比\tenergy\tforce\tvirial\n")
        for info in new_infos:
            record = error_map.get(info.name)
            energy = record.energy if record else float("nan")
            force = record.force if record else float("nan")
            virial = record.virial if record else float("nan")
            handle.write(
                f"{info.name}\t{formatter(info)}\t{energy:.6e}\t{force:.6e}\t{virial:.6e}\n"
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    args = parse_args()
    root = Path.cwd()
    if root.name != "POT":
        print(f"警告：建议在 POT 目录运行，当前目录：{root}")

    iteration_dir = root / args.iteration
    if not iteration_dir.is_dir():
        raise RuntimeError(f"未找到迭代目录：{iteration_dir}")

    dt_dir = iteration_dir / "DT"
    log_step("解析 DT 目录结构...")
    structures, elements = load_structures(dt_dir)

    ss_path = Path(args.ss_path) if args.ss_path else iteration_dir / "ss"
    if not ss_path.is_file():
        raise RuntimeError(f"未找到 ss 文件：{ss_path}")

    log_step("解析 ss 文件...")
    ss_records, thresholds = parse_ss_file(ss_path)

    output_dir = Path(args.output_dir) if args.output_dir else iteration_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    log_step("聚合误差统计...")
    comp_counter, ratio_errors, high_ratio_counts, metric_counts, high_structures = summarize_errors(
        structures,
        ss_records,
        elements,
        thresholds,
    )

    log_step("写出文本报告...")
    write_composition_stats(output_dir / "composition_stats.txt", elements, comp_counter, ratio_errors)
    write_high_error_structures(high_structures, output_dir / "high_error_structures.txt")
    write_sorted_errors(ss_records, "energy", output_dir / "ss_energy.txt")
    write_sorted_errors(ss_records, "force", output_dir / "ss_force.txt")
    write_sorted_errors(ss_records, "virial", output_dir / "ss_virial.txt")

    log_step("生成三角图...")
    figures = []
    fig_all = build_ternary_figure(build_dataset(comp_counter, elements), elements, "全部结构配比分布")
    if fig_all:
        figures.append(fig_all)
    export_html(figures, output_dir / args.html)

    beyond_figs: List[go.Figure] = []
    beyond_figs += build_metric_figures(
        elements,
        "能量 RMSE 超限",
        metric_counts["energy"],
        structures,
        high_structures,
        "over_energy",
    )
    beyond_figs += build_metric_figures(
        elements,
        "受力 RMSE 超限",
        metric_counts["force"],
        structures,
        high_structures,
        "over_force",
    )
    beyond_figs += build_metric_figures(
        elements,
        "位力 RMSE 超限",
        metric_counts["virial"],
        structures,
        high_structures,
        "over_virial",
    )
    export_html(beyond_figs, output_dir / args.high_html)

    log_step("检测新增配比...")
    previous_ratios = collect_previous_ratios(root, elements, args.iteration)
    error_map = {record.name: record for record in ss_records}
    new_infos = [info for info in structures.values() if info.counts not in previous_ratios]
    new_file = output_dir / "new_composition.txt"
    if new_infos:
        write_new_compositions(new_file, new_infos, elements, error_map)
        print(f"新增配比 {len(new_infos)} 条，详情见 {new_file}")
    else:
        new_file.write_text("无新增配比\n", encoding="utf-8")
        print("未发现新增配比。")

    print(f"已生成交互式 HTML：{output_dir / args.html}")
    print(f"高误差配比/结构图：{output_dir / args.high_html}")
    print(f"已输出排序文件：{output_dir / 'ss_energy.txt'}, {output_dir / 'ss_force.txt'}, {output_dir / 'ss_virial.txt'}")
    print(f"配比统计：{output_dir / 'composition_stats.txt'}")
    print(f"高误差结构列表：{output_dir / 'high_error_structures.txt'}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(1)
