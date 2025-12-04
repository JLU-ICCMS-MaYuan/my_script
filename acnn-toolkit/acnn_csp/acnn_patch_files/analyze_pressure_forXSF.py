#!/usr/bin/env python3
"""在 XSF 目录下分析指定迭代的压强结果并输出图表与列表。

示例（在 XSF 目录执行）：
    python analyze_pressure.py IT0 IT1 -p 100 -t 20 --extreme-error 100 --move-xsf 5
    python analyze_pressure.py IT3 -p 80 -t 30  # 仅设置合理阈值
    python analyze_pressure.py -r IT2/extreme_error_structures.txt -m 3  # 直接按文件删除

输出位置：每个迭代的 PNG 与 TXT 文件写入对应的 `IT*` 目录；
--move-xsf 会将极端不合理结构移动到 `XSF/trash_bin`。
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def extract_pressure_from_xsf(xsf_path: Path) -> Optional[float]:
    """从 XSF 文件中解析 VIRIAL 与晶格体积，计算压强 (GPa)。"""

    virial_values: List[float] = []
    lattice: List[List[float]] = []
    try:
        with xsf_path.open("r", encoding="utf-8", errors="ignore") as handle:
            lines = handle.readlines()
    except Exception as exc:  # pragma: no cover - 仅记录异常
        tqdm.write(f"[WARNING] 无法读取 {xsf_path}: {exc}")
        return None

    # 提取 VIRIAL（9 个数字）
    for idx, line in enumerate(lines):
        if line.strip().upper() == "VIRIAL":
            collected: List[float] = []
            j = idx + 1
            while j < len(lines) and len(collected) < 9:
                collected.extend([float(x) for x in lines[j].split()])
                j += 1
            if len(collected) >= 9:
                virial_values = collected[:9]
            break

    # 提取 PRIMVEC（3 行晶格向量）
    for idx, line in enumerate(lines):
        if line.strip().upper() == "PRIMVEC":
            try:
                lattice = [
                    [float(x) for x in lines[idx + 1].split()],
                    [float(x) for x in lines[idx + 2].split()],
                    [float(x) for x in lines[idx + 3].split()],
                ]
            except Exception:  # pragma: no cover - 保守跳过
                lattice = []
            break

    if len(virial_values) < 9 or len(lattice) != 3:
        return None

    a, b, c = (np.array(vec) for vec in lattice)
    volume = abs(np.dot(a, np.cross(b, c)))  # Å^3
    if volume == 0:
        return None

    # virial 单位 eV，压力 = trace(virial)/(3*V)，1 eV/Å^3 = 160.21766208 GPa
    trace_virial = virial_values[0] + virial_values[4] + virial_values[8]
    pressure_gpa = (trace_virial / (3 * volume)) * 160.21766208
    return pressure_gpa


def plot_pressure_distribution(
    pressure_data: List[Dict[str, float]],
    target_pressure_gpa: float,
    threshold_percent: float,
    output_file: Path,
) -> None:
    """绘制主压强分布图并保存。"""

    total_pressures = [d["total_pressure_gpa"] for d in pressure_data]
    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    n, _, _ = ax1.hist(total_pressures, bins=50, alpha=0.7, color="steelblue", edgecolor="black")
    ax1.axvspan(threshold_lower, threshold_upper, alpha=0.15, color="green", label=f"Acceptable (±{threshold_percent}%)", zorder=1)
    ax1.axvline(target_pressure_gpa, color="#E74C3C", linestyle="-", linewidth=3, label=f"Target: {target_pressure_gpa} GPa", zorder=3)
    ax1.axvline(threshold_upper, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Upper: {threshold_upper:.2f} GPa", zorder=2)
    ax1.axvline(threshold_lower, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Lower: {threshold_lower:.2f} GPa", zorder=2)

    y_max = n.max() * 0.95 if len(n) else 1
    ax1.text(target_pressure_gpa, y_max, f"{target_pressure_gpa} GPa", ha="center", va="top", fontsize=11, fontweight="bold", color="#E74C3C")
    ax1.text(threshold_lower, y_max * 0.75, f"{threshold_lower:.2f} GPa", ha="center", va="top", fontsize=9, color="#E67E22")
    ax1.text(threshold_upper, y_max * 0.75, f"{threshold_upper:.2f} GPa", ha="center", va="top", fontsize=9, color="#E67E22")

    ax1.set_xlabel("Total Pressure (GPa)")
    ax1.set_ylabel("Number of Structures")
    ax1.set_title("Pressure Distribution (Histogram)")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    indices = list(range(1, len(total_pressures) + 1))
    colors = ["green" if threshold_lower <= p <= threshold_upper else "red" for p in total_pressures]
    ax2.axhspan(threshold_lower, threshold_upper, alpha=0.15, color="green", label=f"Acceptable (±{threshold_percent}%)", zorder=1)
    ax2.scatter(indices, total_pressures, c=colors, alpha=0.6, s=20, zorder=3)
    ax2.axhline(target_pressure_gpa, color="#E74C3C", linestyle="-", linewidth=3, label=f"Target: {target_pressure_gpa} GPa", zorder=2)
    ax2.axhline(threshold_upper, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Upper: {threshold_upper:.2f} GPa", zorder=2)
    ax2.axhline(threshold_lower, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Lower: {threshold_lower:.2f} GPa", zorder=2)

    ax2.set_xlabel("Structure Index")
    ax2.set_ylabel("Total Pressure (GPa)")
    ax2.set_title("Pressure Distribution (Scatter Plot)")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_extreme_error_distribution(
    extreme_errors: List[Dict[str, float]],
    target_pressure_gpa: float,
    extreme_error_percent: float,
    output_file: Path,
) -> None:
    """绘制极端不合理结构分布图并保存。"""

    extreme_pressures = [e["total_pressure_gpa"] for e in extreme_errors]
    if not extreme_pressures:
        return

    def extract_it_number(path_text: str) -> int:
        match = re.search(r"IT(\d+)", path_text)
        return int(match.group(1)) if match else -1

    for error in extreme_errors:
        it_token_source = str(error.get("outcar_path") or error.get("xsf_path") or "")
        error["it_number"] = extract_it_number(it_token_source)

    extreme_errors_sorted = sorted(extreme_errors, key=lambda x: x["it_number"])
    it_boundaries = []
    current_it = None
    for idx, error in enumerate(extreme_errors_sorted):
        if current_it is not None and error["it_number"] != current_it:
            it_boundaries.append(idx)
        current_it = error["it_number"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
    extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)

    n, _, _ = ax1.hist(extreme_pressures, bins=30, alpha=0.7, color="crimson", edgecolor="black")
    ax1.axvline(target_pressure_gpa, color="#E74C3C", linestyle="-", linewidth=3, label=f"Target: {target_pressure_gpa} GPa", zorder=3)
    ax1.axvline(extreme_upper, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Upper: {extreme_upper:.2f} GPa", zorder=2)
    ax1.axvline(extreme_lower, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Lower: {extreme_lower:.2f} GPa", zorder=2)
    ax1.axvspan(ax1.get_xlim()[0], extreme_lower, alpha=0.15, color="red", zorder=1)
    ax1.axvspan(extreme_upper, ax1.get_xlim()[1], alpha=0.15, color="red", zorder=1)
    y_max = n.max() * 0.95 if len(n) else 1
    ax1.text(target_pressure_gpa, y_max, f"{target_pressure_gpa} GPa", ha="center", va="top", fontsize=11, fontweight="bold", color="#E74C3C")
    ax1.set_xlabel("Total Pressure (GPa)")
    ax1.set_ylabel("Number of Structures")
    ax1.set_title(f"Extreme Error Structures (>{extreme_error_percent}% deviation) - Histogram")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    extreme_pressures_sorted = [e["total_pressure_gpa"] for e in extreme_errors_sorted]
    indices = list(range(1, len(extreme_pressures_sorted) + 1))
    colors = ["darkred" if p < extreme_lower or p > extreme_upper else "orange" for p in extreme_pressures_sorted]
    ax2.axhspan(ax2.get_ylim()[0], extreme_lower, alpha=0.15, color="red", zorder=1)
    ax2.axhspan(extreme_upper, ax2.get_ylim()[1], alpha=0.15, color="red", zorder=1)
    ax2.scatter(indices, extreme_pressures_sorted, c=colors, alpha=0.7, s=30, zorder=3)
    for boundary_idx in it_boundaries:
        ax2.axvline(boundary_idx + 0.5, color="#3498DB", linestyle="--", linewidth=1.5, alpha=0.7, zorder=2)

    boundaries = [0] + it_boundaries + [len(extreme_errors_sorted)]
    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        it_num = extreme_errors_sorted[start].get("it_number", -1)
        if it_num < 0:
            continue
        end = boundaries[i + 1]
        center = (start + end) / 2 + 1
        ax2.text(center, ax2.get_ylim()[1], f"IT{it_num}", ha="center", va="bottom", fontsize=9, color="#3498DB")

    ax2.axhline(target_pressure_gpa, color="#E74C3C", linestyle="-", linewidth=3, label=f"Target: {target_pressure_gpa} GPa", zorder=2)
    ax2.axhline(extreme_upper, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Upper: {extreme_upper:.2f} GPa", zorder=2)
    ax2.axhline(extreme_lower, color="#E67E22", linestyle="--", linewidth=2.5, label=f"Lower: {extreme_lower:.2f} GPa", zorder=2)
    ax2.set_xlabel("Structure Index")
    ax2.set_ylabel("Total Pressure (GPa)")
    ax2.set_title(f"Extreme Error Structures (>{extreme_error_percent}% deviation) - Scatter Plot")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def read_extreme_error_structures_from_file(file_path: Path) -> Optional[List[Dict[str, float]]]:
    """从文件读取极端不合理结构列表。"""

    if not file_path.exists():
        print(f"[ERROR] 文件不存在: {file_path}")
        return None

    extreme_errors: List[Dict[str, float]] = []
    try:
        with file_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                match = re.match(
                    r"\s*\d+\.\s+(.+?)\s+\|\s+([-\d.]+)\s+GPa\s+\|\s+偏差:\s+([-+\d.]+)\s+GPa",
                    line,
                )
                if match:
                    xsf_path = Path(match.group(1))
                    pressure = float(match.group(2))
                    deviation = abs(float(match.group(3)))
                    extreme_errors.append(
                        {
                            "xsf_path": xsf_path,
                            "total_pressure_gpa": pressure,
                            "deviation_gpa": deviation,
                            "file_exists": xsf_path.exists(),
                            "outcar_path": None,
                        }
                    )
        print(f"[INFO] 从 {file_path} 读取到 {len(extreme_errors)} 个极端不合理结构")
        return extreme_errors
    except Exception as exc:
        print(f"[ERROR] 读取 {file_path} 时出错: {exc}")
        return None


def move_extreme_errors_to_trash(
    extreme_errors_to_move: List[Dict[str, float]],
    requested_count: int,
    total_count: int,
    xsf_root: Path,
) -> None:
    """将极端不合理结构移动到 XSF/trash_bin。"""

    print("=" * 70)
    print(f"移动极端不合理结构到垃圾箱 (前 {requested_count} 个，共 {total_count} 个)")
    print("=" * 70)

    trash_dir = xsf_root / "trash_bin"
    trash_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] 垃圾箱目录: {trash_dir.resolve()}")
    print()

    print(f"[INFO] 将移动以下 {len(extreme_errors_to_move)} 个结构（按偏差从大到小）:")
    for idx, error in enumerate(extreme_errors_to_move, 1):
        xsf_path = error.get("xsf_path")
        pressure = error.get("total_pressure_gpa")
        deviation = error.get("deviation_gpa")
        print(f"  {idx:3d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa")
    print()

    moved_count = 0
    not_found_count = 0

    for error in tqdm(extreme_errors_to_move, desc="移动文件", unit="文件"):
        xsf_path = error.get("xsf_path")
        if xsf_path and Path(xsf_path).exists():
            xsf_file = Path(xsf_path)
            it_name = xsf_file.parent.name
            new_name = f"{it_name}_{xsf_file.name}"
            target_path = trash_dir / new_name
            try:
                shutil.move(str(xsf_file), str(target_path))
                moved_count += 1
            except Exception as exc:
                tqdm.write(f"  ✗ 移动失败: {xsf_file}, 错误: {exc}")
        else:
            tqdm.write(f"  ✗ 文件不存在: {xsf_path}")
            not_found_count += 1

    print()
    print(f"[INFO] 成功移动: {moved_count} 个文件")
    if not_found_count > 0:
        print(f"[WARNING] 未找到: {not_found_count} 个文件")
    print("=" * 70)
    print()


def write_outlier_list(
    outliers: List[Dict[str, object]],
    output_file: Path,
    target_pressure_gpa: float,
    threshold_percent: float,
    threshold_lower: float,
    threshold_upper: float,
) -> int:
    outliers.sort(key=lambda x: x["deviation_gpa"], reverse=True)
    missing_outliers_count = 0
    for outlier in outliers:
        xsf_path = outlier.get("xsf_path")
        if xsf_path:
            outlier["file_exists"] = Path(xsf_path).exists()
            if not outlier["file_exists"]:
                missing_outliers_count += 1
        else:
            outlier["file_exists"] = False
            missing_outliers_count += 1

    with output_file.open("w", encoding="utf-8") as handle:
        handle.write(
            f"# 偏离合理范围的结构 (目标压强 {target_pressure_gpa} GPa, 合理阈值 ±{threshold_percent}%)\n"
        )
        handle.write(f"# 共 {len(outliers)} 个结构\n")
        handle.write(f"# 合理范围: {threshold_lower:.2f} ~ {threshold_upper:.2f} GPa\n")
        if missing_outliers_count > 0:
            handle.write(f"# ⚠️ 警告: 有 {missing_outliers_count} 个 XSF 文件不存在（标记为 [文件不存在]）\n")
        handle.write("#\n# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa) | [状态]\n" + "#" + "=" * 68 + "\n\n")

        for idx, outlier in enumerate(outliers, 1):
            xsf_path = outlier.get("xsf_path")
            pressure = outlier.get("total_pressure_gpa")
            deviation = outlier.get("deviation_gpa")
            file_exists = outlier.get("file_exists", False)
            status = "" if file_exists else "  [文件不存在]"
            line = f"{idx:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa{status}\n"
            handle.write(line)
    return missing_outliers_count


def write_extreme_list(
    extreme_errors: List[Dict[str, object]],
    output_file: Path,
    target_pressure_gpa: float,
    extreme_error_percent: float,
) -> int:
    extreme_errors.sort(key=lambda x: x["deviation_gpa"], reverse=True)
    missing_files_count = 0
    for error in extreme_errors:
        xsf_path = error.get("xsf_path")
        if xsf_path:
            error["file_exists"] = Path(xsf_path).exists()
            if not error["file_exists"]:
                missing_files_count += 1
        else:
            error["file_exists"] = False
            missing_files_count += 1

    lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
    upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)

    with output_file.open("w", encoding="utf-8") as handle:
        handle.write(
            f"# 极大误差结构列表 (超过目标压强 {target_pressure_gpa} GPa 的 ±{extreme_error_percent}%)\n"
        )
        handle.write(f"# 共 {len(extreme_errors)} 个结构\n")
        handle.write(f"# 极大误差范围: < {lower:.2f} GPa 或 > {upper:.2f} GPa\n")
        if missing_files_count > 0:
            handle.write(f"# ⚠️ 警告: 有 {missing_files_count} 个 XSF 文件不存在（标记为 [文件不存在]）\n")
        handle.write("#\n# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa) | [状态]\n" + "#" + "=" * 68 + "\n\n")

        for idx, error in enumerate(extreme_errors, 1):
            xsf_path = error.get("xsf_path")
            pressure = error.get("total_pressure_gpa")
            deviation = error.get("deviation_gpa")
            file_exists = error.get("file_exists", False)
            status = "" if file_exists else "  [文件不存在]"
            line = f"{idx:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa{status}\n"
            handle.write(line)
    return missing_files_count


def analyze_iteration(
    it_name: str,
    xsf_root: Path,
    target_pressure_gpa: float,
    threshold_percent: float,
    extreme_error_percent: Optional[float],
    move_extreme_count: Optional[int],
) -> None:
    """针对单个迭代进行压强分析，输出到对应 IT 目录。"""

    it_dir = xsf_root / it_name
    xsf_pattern = it_dir / f"{it_name}_SCF_*"
    xsf_files = [Path(p) for p in glob.glob(str(xsf_pattern)) if Path(p).is_file()]
    if not xsf_files:
        print(f"[ERROR] 未找到 XSF 文件: {xsf_pattern}")
        return

    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    print("=" * 70)
    print(f"迭代 {it_name}：目标压强 {target_pressure_gpa} GPa，合理范围 {threshold_lower:.2f}~{threshold_upper:.2f} GPa")
    if extreme_error_percent is not None:
        extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
        extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)
        print(f"极端不合理范围: < {extreme_lower:.2f} GPa 或 > {extreme_upper:.2f} GPa")
    print(f"共 {len(xsf_files)} 个 XSF，开始提取...")

    pressure_data: List[Dict[str, float]] = []
    outliers: List[Dict[str, float]] = []
    extreme_errors: List[Dict[str, float]] = []

    def handle_pressure(xsf_path: Path, pressure: Optional[float]) -> None:
        if pressure is None:
            return
        entry = {
            "xsf_path": xsf_path,
            "outcar_path": None,
            "total_pressure_gpa": pressure,
            "deviation_gpa": abs(pressure - target_pressure_gpa),
        }

        is_extreme = False
        if extreme_error_percent is not None:
            extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
            extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)
            if pressure < extreme_lower or pressure > extreme_upper:
                is_extreme = True
                extreme_errors.append(entry)

        if not is_extreme:
            pressure_data.append(entry)
            if pressure < threshold_lower or pressure > threshold_upper:
                outliers.append(entry)

    num_workers = max(1, multiprocessing.cpu_count() - 1)
    try:
        if num_workers > 1:
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                future_to_path = {executor.submit(extract_pressure_from_xsf, path): path for path in xsf_files}
                for future in tqdm(as_completed(future_to_path), total=len(xsf_files), desc=f"处理 {it_name}", unit="文件"):
                    xsf_path = future_to_path[future]
                    try:
                        pressure = future.result()
                    except Exception as exc:
                        tqdm.write(f"[WARNING] 处理 {xsf_path} 时出错: {exc}")
                        continue
                    handle_pressure(xsf_path, pressure)
        else:
            for xsf_path in tqdm(xsf_files, desc=f"处理 {it_name}", unit="文件"):
                handle_pressure(xsf_path, extract_pressure_from_xsf(xsf_path))
    except PermissionError:
        tqdm.write("[WARNING] 进程池不可用，改为单进程处理")
        for xsf_path in tqdm(xsf_files, desc=f"处理 {it_name}", unit="文件"):
            handle_pressure(xsf_path, extract_pressure_from_xsf(xsf_path))

    output_dir = xsf_root / it_name
    output_dir.mkdir(parents=True, exist_ok=True)

    if pressure_data:
        total_pressures = [d["total_pressure_gpa"] for d in pressure_data]
        print(f"统计: min={min(total_pressures):.2f}, max={max(total_pressures):.2f}, mean={np.mean(total_pressures):.2f}, median={np.median(total_pressures):.2f}, std={np.std(total_pressures):.2f}")
        plot_pressure_distribution(pressure_data, target_pressure_gpa, threshold_percent, output_dir / "pressure_distribution.png")
    else:
        print(f"[WARNING] {it_name} 未获取到可用于主图的数据")

    if outliers:
        missing_outliers = write_outlier_list(outliers, output_dir / "outlier_structures.txt", target_pressure_gpa, threshold_percent, threshold_lower, threshold_upper)
        print(f"偏离合理范围结构: {len(outliers)}，已写入 {output_dir / 'outlier_structures.txt'}，缺失 {missing_outliers} 个 XSF。")

    if extreme_errors:
        extreme_file = output_dir / "extreme_error_structures.txt"
        missing_extreme = write_extreme_list(extreme_errors, extreme_file, target_pressure_gpa, extreme_error_percent or 0)
        print(f"极大误差结构: {len(extreme_errors)}，已写入 {extreme_file}，缺失 {missing_extreme} 个 XSF。")
        if move_extreme_count:
            to_move = extreme_errors[:move_extreme_count]
            move_extreme_errors_to_trash(to_move, move_extreme_count, len(extreme_errors), xsf_root)
        plot_extreme_error_distribution(extreme_errors, target_pressure_gpa, extreme_error_percent or 0, output_dir / "extreme_error_distribution.png")

    print()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="分析 DFT 结构压强并生成图表（在 XSF 目录运行）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python analyze_pressure.py IT0 IT1 -p 100 -t 20 --extreme-error 100 --move-xsf 5
  python analyze_pressure.py IT3 -p 80 -t 30
  python analyze_pressure.py -r IT2/extreme_error_structures.txt -m 3
        """,
    )

    parser.add_argument("iterations", nargs="*", help="要分析的 IT 目录，如 IT0 IT1")
    parser.add_argument("-p", "--pressure", type=float, default=None, help="目标压强 (GPa) [分析模式必需]")
    parser.add_argument("-t", "--threshold", type=float, default=20.0, help="合理阈值百分比，默认 ±20%%")
    parser.add_argument("-e", "--extreme-error", type=float, default=None, help="极端不合理阈值百分比，超出则单独处理")
    parser.add_argument("-m", "--move-xsf", type=int, default=None, metavar="N", help="移动偏差最大的前 N 个极端不合理结构到 XSF/trash_bin")
    parser.add_argument("-r", "--read-extreme-error-txt", type=str, default=None, metavar="FILE", help="直接读取已有的 extreme_error_structures.txt 并执行移动")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    xsf_root = Path.cwd()
    if xsf_root.name != "XSF":
        print(f"[ERROR] 请在 XSF 目录下运行，当前目录: {xsf_root}")
        return 1

    if args.read_extreme_error_txt:
        if args.move_xsf is None:
            print("[ERROR] 使用 -r 读取模式时必须指定 -m/--move-xsf")
            return 1
        extreme_file = Path(args.read_extreme_error_txt)
        extreme_errors = read_extreme_error_structures_from_file(extreme_file)
        if not extreme_errors:
            print("[ERROR] 未能读取到极端不合理结构")
            return 1
        extreme_errors.sort(key=lambda x: x["deviation_gpa"], reverse=True)
        move_extreme_errors_to_trash(extreme_errors[: args.move_xsf], args.move_xsf, len(extreme_errors), xsf_root)
        print("操作完成！")
        return 0

    if not args.iterations:
        print("[ERROR] 分析模式下需要指定至少一个 IT 目录，例如 IT0 IT1")
        return 1
    if args.pressure is None:
        print("[ERROR] 分析模式下必须指定 -p/--pressure")
        return 1
    if args.move_xsf is not None and args.extreme_error is None:
        print("[ERROR] --move-xsf 需要与 --extreme-error 一起使用")
        return 1

    for it_name in args.iterations:
        analyze_iteration(
            it_name=it_name,
            xsf_root=xsf_root,
            target_pressure_gpa=args.pressure,
            threshold_percent=args.threshold,
            extreme_error_percent=args.extreme_error,
            move_extreme_count=args.move_xsf,
        )

    print("全部处理完成！")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
