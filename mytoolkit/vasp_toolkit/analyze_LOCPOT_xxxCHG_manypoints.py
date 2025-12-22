#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
沿用户指定的多个原子坐标（分数或笛卡尔）依次连线，抽取 LOCPOT / PARCHG / CHGCAR 的 1D 剖面。

用法示例（分数坐标，含周期平移）:
  python analyze_LOCPOT_xxxCHG_manypoints.py \\
    --poscar POSCAR --locpot LOCPOT --parchg PARCHG --chgcar CHGCAR \\
    --direct_points 0.25 0.25 0.25 0 0 0  0 0 0 0 0 0 \\
    --out_prefix line_frac

用法示例（笛卡尔坐标）:
  python analyze_LOCPOT_xxxCHG_manypoints.py \\
    --cart_points C 0.25 0.25 0.25  C 0.0 0.0 0.0 \\
    --out_prefix line_cart
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# -----------------------------
# Parsing POSCAR-like header
# -----------------------------
def read_poscar_header(lines):
    idx = 0
    _ = lines[idx].rstrip("\n"); idx += 1
    scale = float(lines[idx].split()[0]); idx += 1

    lattice = []
    for _ in range(3):
        lattice.append([float(x) for x in lines[idx].split()[:3]])
        idx += 1
    lattice = np.array(lattice, dtype=float) * scale  # Å

    tokens = lines[idx].split()
    def _all_int(tok_list):
        try:
            _ = [int(t) for t in tok_list]
            return True
        except Exception:
            return False

    if _all_int(tokens):
        symbols = None
        counts = [int(t) for t in tokens]
        idx += 1
    else:
        symbols = tokens
        idx += 1
        counts = [int(t) for t in lines[idx].split()]
        idx += 1

    line = lines[idx].strip().lower()
    if line.startswith("s"):
        idx += 1
        line = lines[idx].strip().lower()
    coord_type = line
    idx += 1

    n_atoms = int(np.sum(counts))
    positions = []
    for _ in range(n_atoms):
        parts = lines[idx].split()
        positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
        idx += 1

    while idx < len(lines) and len(lines[idx].strip()) == 0:
        idx += 1

    return lattice, symbols, counts, coord_type, np.array(positions, dtype=float), idx


def expand_atom_symbols(symbols, counts):
    if symbols is None:
        symbols = [f"El{i + 1}" for i in range(len(counts))]
    expanded = []
    for sym, cnt in zip(symbols, counts):
        expanded.extend([sym] * cnt)
    return expanded


def get_frac_positions_from_poscar(lines, lattice):
    _, symbols, counts, coord_type, positions, _ = read_poscar_header(lines)
    coord = coord_type.strip().lower()
    if coord.startswith("d"):
        frac_positions = positions
    elif coord.startswith("c") or coord.startswith("k"):
        frac_positions = np.array([cart_to_frac(lattice, p) for p in positions])
    else:
        raise ValueError(f"Unknown coordinate type: {coord_type}")
    atom_symbols = expand_atom_symbols(symbols, counts)
    return frac_positions, atom_symbols


def read_vasp_volumetric(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    lattice, symbols, counts, coord_type, pos, idx = read_poscar_header(lines)
    nx, ny, nz = [int(x) for x in lines[idx].split()[:3]]
    idx += 1
    ngrid = nx * ny * nz
    vals = []
    while idx < len(lines) and len(vals) < ngrid:
        parts = lines[idx].split()
        if parts:
            vals.extend([float(x) for x in parts])
        idx += 1
    if len(vals) < ngrid:
        raise RuntimeError(f"{filepath}: not enough grid values. got {len(vals)}, need {ngrid}")
    data = np.array(vals[:ngrid], dtype=float).reshape((nx, ny, nz), order="F")
    return lattice, (nx, ny, nz), data


# -----------------------------
# Geometry helpers
# -----------------------------
def cart_to_frac(lattice, r_cart):
    A = lattice.T
    return np.linalg.solve(A, r_cart)


def frac_to_cart(lattice, f):
    A = lattice.T
    return A @ f


def sample_trilinear_periodic(data, frac, grid):
    nx, ny, nz = grid
    f = frac - np.floor(frac)
    x = f[0] * nx
    y = f[1] * ny
    z = f[2] * nz
    i0 = int(np.floor(x)) % nx
    j0 = int(np.floor(y)) % ny
    k0 = int(np.floor(z)) % nz
    i1 = (i0 + 1) % nx
    j1 = (j0 + 1) % ny
    k1 = (k0 + 1) % nz
    tx = x - np.floor(x)
    ty = y - np.floor(y)
    tz = z - np.floor(z)
    c000 = data[i0, j0, k0]; c100 = data[i1, j0, k0]
    c010 = data[i0, j1, k0]; c110 = data[i1, j1, k0]
    c001 = data[i0, j0, k1]; c101 = data[i1, j0, k1]
    c011 = data[i0, j1, k1]; c111 = data[i1, j1, k1]
    c00 = c000 * (1 - tx) + c100 * tx
    c10 = c010 * (1 - tx) + c110 * tx
    c01 = c001 * (1 - tx) + c101 * tx
    c11 = c011 * (1 - tx) + c111 * tx
    c0 = c00 * (1 - ty) + c10 * ty
    c1 = c01 * (1 - ty) + c11 * ty
    return c0 * (1 - tz) + c1 * tz


def make_cell_edges():
    corners = np.array([[i, j, k] for i in [0, 1] for j in [0, 1] for k in [0, 1]], dtype=float)
    edges = []
    for a in range(len(corners)):
        for b in range(a + 1, len(corners)):
            diff = np.abs(corners[a] - corners[b])
            if np.isclose(np.sum(diff), 1.0) and np.count_nonzero(diff) == 1:
                edges.append((corners[a], corners[b]))
    return edges


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max([x_range, y_range, z_range])
    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)
    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2])


def plot_polyline_geometry(lattice, points_frac, out_png):
    points_cart = np.array([frac_to_cart(lattice, p) for p in points_frac])
    edges = make_cell_edges()
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    for f0, f1 in edges:
        p0 = frac_to_cart(lattice, f0)
        p1 = frac_to_cart(lattice, f1)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], color="0.6", linewidth=1)
    ax.plot(points_cart[:, 0], points_cart[:, 1], points_cart[:, 2], color="tab:red", linewidth=2, label="path")
    ax.scatter(points_cart[:, 0], points_cart[:, 1], points_cart[:, 2], color="tab:green", s=30)
    for idx, p in enumerate(points_cart, start=1):
        ax.text(p[0], p[1], p[2], f"P{idx}", color="tab:green", fontsize=8)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.legend()
    set_axes_equal(ax)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# -----------------------------
# Sampling along polyline
# -----------------------------
def compute_grid_step(lattice, grid):
    nx, ny, nz = grid
    a1, a2, a3 = lattice
    return min(np.linalg.norm(a1) / nx, np.linalg.norm(a2) / ny, np.linalg.norm(a3) / nz)


def resample_polyline(points_cart, step=None, npts=None):
    seg_vecs = points_cart[1:] - points_cart[:-1]
    seg_len = np.linalg.norm(seg_vecs, axis=1)
    total_len = float(np.sum(seg_len))
    if total_len < 1e-12:
        raise ValueError("Total polyline length is zero.")
    if npts is not None:
        if npts < 2:
            raise ValueError("npts must be >= 2.")
        step = total_len / (npts - 1)
    else:
        if step is None or step <= 0:
            raise ValueError("Invalid step.")
        npts = int(np.floor(total_len / step)) + 1
        if npts < 2:
            npts = 2
    s_target = np.linspace(0.0, total_len, npts)
    cum_len = np.concatenate([[0.0], np.cumsum(seg_len)])
    pts = np.zeros((npts, 3), dtype=float)
    j = 0
    for i, st in enumerate(s_target):
        while j < len(seg_len) - 1 and st > cum_len[j + 1] + 1e-12:
            j += 1
        if seg_len[j] < 1e-12:
            pts[i] = points_cart[j]
            continue
        t = (st - cum_len[j]) / seg_len[j]
        pts[i] = points_cart[j] + t * seg_vecs[j]
    return s_target, pts


def extract_profile_along_polyline(lattice, grid, V, rho_parchg, rho_chg, points_frac, step, npts):
    points_cart = np.array([frac_to_cart(lattice, p) for p in points_frac], dtype=float)
    if step is None and npts is None:
        step = compute_grid_step(lattice, grid)
    s, pts_cart = resample_polyline(points_cart, step=step, npts=npts)
    V_line = np.zeros_like(s)
    rho_line = np.zeros_like(s)
    rho_chg_line = np.zeros_like(s) if rho_chg is not None else None
    for idx, rc in enumerate(pts_cart):
        f = cart_to_frac(lattice, rc)
        V_line[idx] = sample_trilinear_periodic(V, f, grid)
        rho_line[idx] = sample_trilinear_periodic(rho_parchg, f, grid)
        if rho_chg_line is not None:
            rho_chg_line[idx] = sample_trilinear_periodic(rho_chg, f, grid)
    return s, V_line, rho_line, rho_chg_line


def compute_vertex_s(points_cart):
    """返回各 polyline 顶点在总长坐标中的 s 值（Å），用于标注。"""
    seg_vecs = points_cart[1:] - points_cart[:-1]
    seg_len = np.linalg.norm(seg_vecs, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    return cum


# -----------------------------
# Input parsing for points
# -----------------------------
def parse_cart_points(cart_args, lattice):
    if len(cart_args) % 4 != 0:
        raise ValueError("Cart points must be groups of: Elem x y z")
    points_frac = []
    labels = []
    for i in range(0, len(cart_args), 4):
        elem = cart_args[i]
        x, y, z = map(float, cart_args[i + 1:i + 4])
        frac = cart_to_frac(lattice, np.array([x, y, z], dtype=float))
        points_frac.append(frac)
        labels.append(elem)
    return np.array(points_frac, dtype=float), labels


def parse_direct_points(direct_args):
    if len(direct_args) % 7 != 0:
        raise ValueError("Direct points must be groups of: Elem fx fy fz sx sy sz")
    points_frac = []
    labels = []
    for i in range(0, len(direct_args), 7):
        elem = direct_args[i]
        fx, fy, fz = map(float, direct_args[i + 1:i + 4])
        sx, sy, sz = map(float, direct_args[i + 4:i + 7])
        frac = np.array([fx + sx, fy + sy, fz + sz], dtype=float)
        points_frac.append(frac)
        labels.append(str(elem))
    return np.array(points_frac, dtype=float), labels


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poscar", default="POSCAR", help="POSCAR file")
    ap.add_argument("--locpot", default="LOCPOT", help="LOCPOT file")
    ap.add_argument("--parchg", default="PARCHG", help="PARCHG file")
    ap.add_argument("--chgcar", default="CHGCAR", help="CHGCAR file (default: CHGCAR, auto-skip if missing)")
    ap.add_argument("--cart_points", "-c", nargs="+", help="Cartesian points as: Elem x y z (repeatable)")
    ap.add_argument("--direct_points", "-d", nargs="+", help="Direct points as: Elem fx fy fz sx sy sz (repeatable)")
    ap.add_argument("--step", type=float, default=None, help="Sampling step (Å) along the polyline")
    ap.add_argument("--npts", type=int, default=None, help="Number of sampling points (overrides step)")
    ap.add_argument("--out_prefix", default="line_profile", help="Output prefix")
    args = ap.parse_args()

    if args.cart_points is None and args.direct_points is None:
        raise ValueError("Please provide --cart_points or --direct_points.")

    with open(args.poscar, "r") as f:
        poscar_lines = f.readlines()
    lattice_poscar, *_ = read_poscar_header(poscar_lines)
    lattice = lattice_poscar

    lattice_L, grid_L, V = read_vasp_volumetric(args.locpot)
    lattice_C, grid_C, rho = read_vasp_volumetric(args.parchg)
    rho_chg = None
    lattice_G = grid_G = None
    if args.chgcar and os.path.exists(args.chgcar):
        lattice_G, grid_G, rho_chg = read_vasp_volumetric(args.chgcar)
    elif args.chgcar:
        print(f"CHGCAR not found at {args.chgcar}, skip plotting CHGCAR.")

    if grid_L != grid_C:
        raise RuntimeError(f"Grid mismatch: LOCPOT {grid_L} vs PARCHG {grid_C}")
    if rho_chg is not None and grid_L != grid_G:
        raise RuntimeError(f"Grid mismatch: LOCPOT {grid_L} vs CHGCAR {grid_G}")
    grid = grid_L

    cell_volume = float(abs(np.linalg.det(lattice)))
    rho = rho / cell_volume
    rho_unit = "Å^-3"
    rho_integral = float(rho.mean() * cell_volume)
    print(f"PARCHG converted to density ({rho_unit}); integral ≈ {rho_integral:.10g} e")
    rho_chg_line = None
    if rho_chg is not None:
        rho_chg = rho_chg / cell_volume
        rho_chg_integral = float(rho_chg.mean() * cell_volume)
        print(f"CHGCAR converted to density ({rho_unit}); integral ≈ {rho_chg_integral:.10g} e")

    if args.cart_points is not None:
        points_frac, labels = parse_cart_points(args.cart_points, lattice)
    else:
        points_frac, labels = parse_direct_points(args.direct_points)
    if len(points_frac) < 2:
        raise ValueError("At least two points are required to define a polyline.")

    s, V_line, rho_line, rho_chg_line = extract_profile_along_polyline(
        lattice=lattice, grid=grid, V=V, rho_parchg=rho, rho_chg=rho_chg,
        points_frac=points_frac, step=args.step, npts=args.npts
    )

    # 顶点在路径上的 s 位置，用于标注
    points_cart_for_mark = np.array([frac_to_cart(lattice, p) for p in points_frac])
    s_vertices = compute_vertex_s(points_cart_for_mark)

    print(f"LOCPOT along path: V_min = {V_line.min():.6f} eV, V_max = {V_line.max():.6f} eV")

    geom_png = f"{args.out_prefix}_geom.png"
    try:
        plot_polyline_geometry(lattice, points_frac, geom_png)
        print("Saved geometry figure:", geom_png)
    except Exception as e:
        print("Geometry plot failed:", e)

    out_dat = f"{args.out_prefix}.dat"
    if rho_chg_line is not None:
        header = f"s(Angstrom)\tV_locpot(eV)\trho_parchg({rho_unit})\trho_chgcar({rho_unit})"
        data_out = np.column_stack([s, V_line, rho_line, rho_chg_line])
    else:
        header = f"s(Angstrom)\tV_locpot(eV)\trho_parchg({rho_unit})"
        data_out = np.column_stack([s, V_line, rho_line])
    np.savetxt(out_dat, data_out, header=header, delimiter="\t", fmt="%.10e")

    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    color_v = "tab:blue"
    color_rho = "tab:orange"
    ax1.plot(s, V_line, label="V (LOCPOT)", color=color_v)
    ax1.set_xlabel("s (Å) along path")
    ax1.set_ylabel("V (eV)", color=color_v)
    ax1.tick_params(axis="y", labelcolor=color_v)
    ax1.grid(True)
    ax2 = ax1.twinx()
    ax2.plot(s, rho_line, label="rho (PARCHG)", color=color_rho)
    if rho_chg_line is not None:
        ax2.plot(s, rho_chg_line, label="rho (CHGCAR)", color="tab:green")
    ax2.set_ylabel(f"Density ({rho_unit})", color=color_rho)
    ax2.tick_params(axis="y", labelcolor=color_rho)

    # 标注顶点位置
    label_dx = compute_grid_step(lattice, grid)
    for s_i, lab in zip(s_vertices, labels):
        ax1.axvline(s_i, color="0.5", linestyle="--", linewidth=1, alpha=0.7)
        ax1.text(s_i, 1.02, lab, transform=ax1.get_xaxis_transform(),
                 ha="center", va="bottom", color="0.2", fontsize=8, fontweight="bold")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="best")

    out_png = f"{args.out_prefix}.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print("Done.")
    print("Saved:", out_dat)
    print("Saved:", out_png)


if __name__ == "__main__":
    main()
