#!/usr/bin/env python3
"""Python equivalent of gam.f90 to generate bands and spoint outputs."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

RY_TO_CM1 = 109_737.57990
RY_TO_THZ = 3.289_828e3
PI = math.pi
DEFAULT_SPOINT_INTERVAL = 21


def parse_header(line: str) -> Tuple[int, int]:
    """Extract nband and nks from the first line of the frequency file."""
    nband_match = re.search(r"nbnd=\s*(\d+)", line, flags=re.IGNORECASE)
    nks_match = re.search(r"nks=\s*(\d+)", line, flags=re.IGNORECASE)
    if not nband_match or not nks_match:
        raise ValueError("nbnd or nks not found in the frequency file header")
    return int(nband_match.group(1)), int(nks_match.group(1))


def _read_values(
    lines: Sequence[str], start_idx: int, expected: int
) -> Tuple[List[float], int]:
    """Read floating numbers across lines until reaching expected count."""
    values: List[float] = []
    idx = start_idx
    while len(values) < expected and idx < len(lines):
        chunk = lines[idx].strip()
        idx += 1
        if not chunk:
            continue
        values.extend(float(item) for item in chunk.split())
    if len(values) < expected:
        raise ValueError("not enough frequency values; file may be incomplete")
    return values[:expected], idx


def read_frequency_file(
    path: Path,
) -> Tuple[int, int, List[Tuple[float, float, float]], List[List[float]]]:
    """Load k-points and frequencies from the QE .freq file."""
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError("frequency file is empty")

    nband, nks = parse_header(lines[0])
    idx = 1
    kpoints: List[Tuple[float, float, float]] = []
    frequencies: List[List[float]] = []

    for _ in range(nks):
        # Read k-point coordinates.
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        if idx >= len(lines):
            raise ValueError("frequency file ended while reading k-points")
        k_parts = lines[idx].split()
        if len(k_parts) != 3:
            raise ValueError(f"expected 3 k-point components, got: {lines[idx]!r}")
        kpoints.append(tuple(float(item) for item in k_parts))
        idx += 1

        # Read frequency values for this k-point.
        values, idx = _read_values(lines, idx, nband)
        frequencies.append(values)

    return nband, nks, kpoints, frequencies


def read_gamma_from_gam_lines(
    lines: Sequence[str], nband: int, nks: int, broadening: float | None
) -> List[List[float]]:
    """Parse gam.lines format with multiple broadening blocks."""
    markers: List[Tuple[int, float]] = []
    for idx, line in enumerate(lines):
        m = re.search(r"Broadening\s+([0-9.]+)", line, re.IGNORECASE)
        if m:
            markers.append((idx, float(m.group(1))))

    if not markers:
        raise ValueError("gam.lines 未找到 Broadening 段")

    preferred = broadening if broadening is not None else 0.02
    target_idx = None
    for idx, value in markers:
        if preferred is not None and abs(value - preferred) < 1e-9:
            target_idx = idx
            break
    if target_idx is None:
        raise ValueError(
            f"未找到指定 Broadening={preferred} 的数据块；可使用 --broadening 选择现有值"
        )

    # collect lines after target_idx until next Broadening or EOF
    gammas: List[List[float]] = []
    for line in lines[target_idx + 1 :]:
        if re.search(r"Broadening", line, re.IGNORECASE):
            break
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        values = [float(item) for item in parts[1:]]
        if len(values) < nband:
            raise ValueError(f"gamma row has {len(values)} values, expected {nband}")
        gammas.append(values[:nband])
        if len(gammas) >= nks:
            break

    if len(gammas) != nks:
        raise ValueError(
            f"gamma rows {len(gammas)} shorter than expected {nks} for selected broadening"
        )
    return gammas


def read_gamma_file(
    path: Path, nband: int, nks: int, broadening: float | None
) -> List[List[float]]:
    """Load gamma from gam.lines (Broadening blocks)."""
    lines = path.read_text().splitlines()
    return read_gamma_from_gam_lines(lines, nband, nks, broadening)


def accumulate_kpath(kpoints: Sequence[Sequence[float]]) -> List[float]:
    """Compute cumulative distance along the k-path."""
    totals: List[float] = []
    total = 0.0
    previous = None
    for coords in kpoints:
        if previous is not None:
            dx = sum((c - p) ** 2 for c, p in zip(coords, previous))
            total += math.sqrt(dx)
        totals.append(total)
        previous = coords
    return totals


def compute_lambda_and_gamma(
    frequencies: Sequence[Sequence[float]],
    gammas: Sequence[Sequence[float]],
    dosfit: float,
) -> Tuple[List[List[float]], List[List[float]]]:
    """Compute lambda and adjusted gamma arrays."""
    lambda_rows: List[List[float]] = []
    gamma_rows: List[List[float]] = []
    for freq_row, gamma_row in zip(frequencies, gammas):
        lambda_band: List[float] = []
        gamma_band: List[float] = []
        for freq, gamma in zip(freq_row, gamma_row):
            if gamma > 0:
                lamb = gamma / PI / RY_TO_THZ / dosfit / (freq / RY_TO_CM1) ** 2
                lamb = math.sqrt(lamb / PI)
                lambda_band.append(lamb)
                gamma_band.append(gamma)
            else:
                lambda_band.append(0.0)
                gamma_band.append(abs(gamma))
        lambda_rows.append(lambda_band)
        gamma_rows.append(gamma_band)
    return lambda_rows, gamma_rows


def write_bands(
    path: Path,
    x_totals: Sequence[float],
    frequencies: Sequence[Sequence[float]],
    lambdas: Sequence[Sequence[float]],
    gammas: Sequence[Sequence[float]],
) -> None:
    """Write bands.dat-style output."""
    nband = len(frequencies[0])
    with path.open("w") as handle:
        for band_idx in range(nband):
            for x, freq_row, lambda_row, gamma_row in zip(
                x_totals, frequencies, lambdas, gammas
            ):
                handle.write(
                    f"{x:12.6f}  {freq_row[band_idx]:12.6f}  "
                    f"{lambda_row[band_idx]:12.6f}  {gamma_row[band_idx]:12.6f}\n"
                )
            handle.write("\n")


def collect_spoints(x_totals: Sequence[float], interval: int) -> List[float]:
    """Collect cumulative distances at high-symmetry boundaries."""
    points: List[float] = [x_totals[0]]
    for idx, x in enumerate(x_totals, start=1):
        if idx % interval == 0:
            points.append(x)
    if points[-1] != x_totals[-1]:
        points.append(x_totals[-1])
    return points


def write_spoints(
    path: Path, points: Sequence[float], labels: Sequence[str]
) -> None:
    """Write labeled high-symmetry points: label  distance."""
    with path.open("w") as handle:
        for label, x in zip(labels, points):
            handle.write(f"{label}\t{x:10.4f}\n")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate bands.dat and spoint.dat (Python replacement of gam.f90)",
    )
    parser.add_argument(
        "-f",
        "--freq-file",
        default="SH3.freq",
        type=Path,
        help="frequency file produced by matdyn.in_freq",
    )
    parser.add_argument(
        "-g",
        "--gamma-file",
        default="gam.lines",
        type=Path,
        help="gamma 数据文件，要求 gam.lines 格式（包含 Broadening 块）",
    )
    parser.add_argument(
        "-b",
        "--broadening",
        type=float,
        default=None,
        help="choose Broadening value in gam.lines (default 0.02 if available)",
    )
    parser.add_argument(
        "-d",
        "--dosfit",
        required=True,
        type=float,
        help="dosfit at (0,0,0), e.g. 3.497372",
    )
    parser.add_argument(
        "--bands-output",
        default="bands.dat",
        type=Path,
        help="bands output path, default bands.dat",
    )
    parser.add_argument(
        "--spoint-output",
        default="spoint.dat",
        type=Path,
        help="spoint output path, default spoint.dat",
    )
    parser.add_argument(
        "--spoint-interval",
        default=DEFAULT_SPOINT_INTERVAL,
        type=int,
        help="write spoint every N k-points (Fortran default is 21)",
    )
    parser.add_argument(
        "--spoint-labels",
        type=str,
        default=None,
        help="comma-separated labels for high-symmetry points (must match point count)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    if args.dosfit <= 0:
        raise ValueError("dosfit 必须为正数")

    nband, nks, kpoints, frequencies = read_frequency_file(args.freq_file)
    gammas = read_gamma_file(args.gamma_file, nband, nks, args.broadening)
    x_totals = accumulate_kpath(kpoints)
    lambdas, adjusted_gammas = compute_lambda_and_gamma(
        frequencies, gammas, args.dosfit
    )

    write_bands(args.bands_output, x_totals, frequencies, lambdas, adjusted_gammas)

    spoints = collect_spoints(x_totals, args.spoint_interval)
    if args.spoint_labels:
        labels = [item.strip() for item in args.spoint_labels.split(",") if item.strip()]
        if len(labels) != len(spoints):
            raise ValueError(
                f"提供的高对称点标签数量({len(labels)})与点数({len(spoints)})不一致"
            )
    else:
        labels = [f"P{i+1}" for i in range(len(spoints))]
    write_spoints(args.spoint_output, spoints, labels)


if __name__ == "__main__":
    main()
