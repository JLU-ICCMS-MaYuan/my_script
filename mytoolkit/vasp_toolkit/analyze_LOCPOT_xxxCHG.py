#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract 1D profile of LOCPOT, PARCHG (and optionally CHGCAR) along a specified fractional direction.

Usage example:
  python analyze_LOCPOT_xxxCHG.py \
    --poscar POSCAR \
    --locpot LOCPOT \
    --parchg PARCHG \
    --direction 1 0 0 \\    # fractional direction [u v w]
    --start 0.0 0.0 0.0 \\  # fractional start point
    --length 8.0 \\         # optional, default: to boundary
    --out_prefix line_1_0_0_from_0_0_0

Outputs:
  <out_prefix>.dat   (s[Å], V(eV), rho)
  <out_prefix>.png

Line definition:
  use --direction + --start (fractional)

Atom markers:
  Atoms within a tolerance to the line are marked with vertical dashed lines and labels.
  Length beyond the unit cell is allowed; sampling uses periodic boundary conditions.
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
    """
    Parse POSCAR/CHGCAR/LOCPOT-like header.
    Return:
      lattice (3x3) in Å, element symbols list, counts list, coord_type, positions (Nx3 fractional), next_line_index
    """
    idx = 0
    comment = lines[idx].rstrip("\n"); idx += 1
    scale = float(lines[idx].split()[0]); idx += 1

    lattice = []
    for _ in range(3):
        lattice.append([float(x) for x in lines[idx].split()[:3]])
        idx += 1
    lattice = np.array(lattice, dtype=float) * scale  # Å

    # Next line can be element symbols OR counts (VASP4 vs VASP5)
    tokens = lines[idx].split()
    # if all tokens are ints => VASP4, else VASP5
    def _all_int(tok_list):
        try:
            _ = [int(t) for t in tok_list]
            return True
        except Exception:
            return False

    if _all_int(tokens):
        # VASP4: no symbols line
        symbols = None
        counts = [int(t) for t in tokens]
        idx += 1
    else:
        symbols = tokens
        idx += 1
        counts = [int(t) for t in lines[idx].split()]
        idx += 1

    # Optional "Selective dynamics"
    line = lines[idx].strip().lower()
    selective = False
    if line.startswith("s"):
        selective = True
        idx += 1
        line = lines[idx].strip().lower()

    coord_type = line  # direct/cartesian
    idx += 1

    n_atoms = int(np.sum(counts))
    positions = []
    for _ in range(n_atoms):
        parts = lines[idx].split()
        positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
        idx += 1
    positions = np.array(positions, dtype=float)

    # skip blank lines
    while idx < len(lines) and len(lines[idx].strip()) == 0:
        idx += 1

    return lattice, symbols, counts, coord_type, positions, idx


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
    positions = np.array(positions, dtype=float)

    if coord.startswith("d"):
        frac_positions = positions
    elif coord.startswith("c") or coord.startswith("k"):
        frac_positions = np.array([cart_to_frac(lattice, p) for p in positions])
    else:
        raise ValueError(f"Unknown coordinate type: {coord_type}")

    atom_symbols = expand_atom_symbols(symbols, counts)
    return frac_positions, atom_symbols


def read_vasp_volumetric(filepath):
    """
    Read VASP volumetric file: LOCPOT / CHGCAR / PARCHG (same structure up to first dataset).
    Return:
      lattice(3x3 Å), grid (Nx,Ny,Nz), data (Nx,Ny,Nz) as float
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    lattice, symbols, counts, coord_type, pos, idx = read_poscar_header(lines)

    # Grid line: Nx Ny Nz
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

    vals = np.array(vals[:ngrid], dtype=float)

    # VASP writes x fastest, then y, then z => reshape with Fortran order
    data = vals.reshape((nx, ny, nz), order="F")

    return lattice, (nx, ny, nz), data


# -----------------------------
# Geometry helpers
# -----------------------------
def cart_to_frac(lattice, r_cart):
    """Convert Cartesian (Å) to fractional."""
    A = lattice.T  # columns a1,a2,a3
    f = np.linalg.solve(A, r_cart)
    return f


def frac_to_cart(lattice, f):
    """Convert fractional to Cartesian (Å)."""
    A = lattice.T
    return A @ f


def sample_trilinear_periodic(data, frac, grid):
    """Trilinear interpolation with periodic boundary."""
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

    c000 = data[i0, j0, k0]
    c100 = data[i1, j0, k0]
    c010 = data[i0, j1, k0]
    c110 = data[i1, j1, k0]
    c001 = data[i0, j0, k1]
    c101 = data[i1, j0, k1]
    c011 = data[i0, j1, k1]
    c111 = data[i1, j1, k1]

    c00 = c000 * (1 - tx) + c100 * tx
    c10 = c010 * (1 - tx) + c110 * tx
    c01 = c001 * (1 - tx) + c101 * tx
    c11 = c011 * (1 - tx) + c111 * tx

    c0 = c00 * (1 - ty) + c10 * ty
    c1 = c01 * (1 - ty) + c11 * ty

    c = c0 * (1 - tz) + c1 * tz
    return c


# -----------------------------
# Visualization helpers
# -----------------------------
def make_cell_edges():
    """Return list of edge endpoints (in fractional coords) for the unit cell."""
    corners = np.array([[i, j, k] for i in [0, 1] for j in [0, 1] for k in [0, 1]], dtype=float)
    edges = []
    for a in range(len(corners)):
        for b in range(a + 1, len(corners)):
            diff = np.abs(corners[a] - corners[b])
            if np.isclose(np.sum(diff), 1.0) and np.count_nonzero(diff) == 1:
                edges.append((corners[a], corners[b]))
    return edges


def set_axes_equal(ax):
    """Equalize 3D axes limits for better geometry view."""
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


def plot_line_geometry(lattice, start_frac, direction_frac, length, out_png):
    """Draw unit cell, start point, and the line segment in Cartesian space."""
    start_frac = np.array(start_frac, dtype=float)
    direction_frac = np.array(direction_frac, dtype=float)

    d_cart = frac_to_cart(lattice, direction_frac)
    d_norm = np.linalg.norm(d_cart)
    if d_norm < 1e-12:
        raise ValueError("Direction vector is zero in Cartesian space.")

    start_cart = frac_to_cart(lattice, start_frac)
    end_cart = start_cart + d_cart / d_norm * length

    edges = make_cell_edges()
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")

    # unit cell edges
    for f0, f1 in edges:
        p0 = frac_to_cart(lattice, f0)
        p1 = frac_to_cart(lattice, f1)
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], color="0.6", linewidth=1)

    # start and line segment
    ax.scatter([start_cart[0]], [start_cart[1]], [start_cart[2]], color="tab:green", s=40, label="start")
    ax.plot([start_cart[0], end_cart[0]], [start_cart[1], end_cart[1]], [start_cart[2], end_cart[2]],
            color="tab:red", linewidth=2, label="line")
    ax.text(end_cart[0], end_cart[1], end_cart[2], "direction", color="tab:red", fontsize=8)

    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.legend()
    set_axes_equal(ax)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# -----------------------------
# Main: line profile along given direction
# -----------------------------
def compute_line_max_length(lattice, start_frac, direction_frac):
    """Return distance (Å) to the nearest unit-cell boundary along the given direction (periodic-safe, sign-insensitive)."""
    start_frac = np.array(start_frac, dtype=float)
    direction_frac = np.array(direction_frac, dtype=float)

    if np.linalg.norm(direction_frac) < 1e-12:
        raise ValueError("Direction vector is zero.")

    # ensure inside [0,1) for boundary distance calculation
    start_frac = start_frac - np.floor(start_frac)

    t_candidates = []
    for i in range(3):
        di = direction_frac[i]
        if abs(di) < 1e-12:
            continue
        # distances to both faces (0 and 1) in this dimension; take positive intersections
        for face in (0.0, 1.0):
            delta = abs(face - start_frac[i])
            if delta > 1e-12:
                t_candidates.append(delta / abs(di))

    if not t_candidates:
        raise ValueError("Direction does not intersect unit-cell boundaries.")

    t_hit = min(t_candidates)  # nearest boundary along given direction (sign-aware, periodic)
    d_cart = frac_to_cart(lattice, direction_frac)
    return np.linalg.norm(d_cart) * t_hit


def compute_grid_step(lattice, grid):
    """Estimate grid spacing (Å) using VASP grid."""
    nx, ny, nz = grid
    a1, a2, a3 = lattice
    steps = [
        np.linalg.norm(a1) / nx,
        np.linalg.norm(a2) / ny,
        np.linalg.norm(a3) / nz,
    ]
    return min(steps)


def extract_profile_along_direction(lattice, grid, V, rho, direction_frac, start_frac, length, step, npts, rho_chg=None):
    """
    Extract profile along a specified fractional direction [u, v, w] from a fractional start point.
    """
    direction_frac = np.array(direction_frac, dtype=float)  # 将方向向量转为 float 数组
    start_frac = np.array(start_frac, dtype=float)          # 将起点转为 float 数组（分数坐标）

    d_cart = frac_to_cart(lattice, direction_frac)          # 方向向量转换到笛卡尔坐标
    d_norm = np.linalg.norm(d_cart)                         # 笛卡尔方向的模长
    if d_norm < 1e-12:                                      # 模长过小视为零向量
        raise ValueError("Direction vector is zero in Cartesian space.")

    max_length = compute_line_max_length(lattice, start_frac, direction_frac)  # 计算到晶胞边界的距离

    if length is None:                                      # 若用户未指定长度
        length = max_length                                 # 使用到边界的距离（方向符号已考虑）
    else:
        if length <= 0:                                     # 用户指定长度必须为正
            raise ValueError("Length must be positive.")
        if length > max_length + 1e-8:                      # 超出边界则提示周期采样
            print(f"Warning: length {length:.6f} Å exceeds boundary {max_length:.6f} Å, using periodic sampling.")

    if npts is not None:                                    # 如果用户指定点数
        if npts < 2:
            raise ValueError("npts must be >= 2.")
        step = length / (npts - 1)                          # 用点数反算步长
    else:
        if step is None:                                    # 用户未指定步长
            step = compute_grid_step(lattice, grid)         # 默认用最小网格步长
        if step <= 0:
            raise ValueError("Step must be positive.")
        npts = int(np.floor(length / step)) + 1             # 由步长算点数
        if npts < 2:                                        # 至少两个点
            npts = 2

    s = np.linspace(0.0, length, npts)                      # 沿线距离数组（单位 Å）
    V_line = np.zeros_like(s)                               # 存放 LOCPOT 剖面
    rho_line = np.zeros_like(s)                             # 存放 PARCHG 剖面
    rho_chg_line = np.zeros_like(s) if rho_chg is not None else None  # 存放 CHGCAR 剖面（可选）

    r0 = frac_to_cart(lattice, start_frac)                  # 起点转为笛卡尔坐标
    d_unit = d_cart / d_norm                                # 单位方向向量

    for idx, si in enumerate(s):                            # 遍历每个采样距离
        r = r0 + si * d_unit                                # 沿线的笛卡尔坐标
        f = cart_to_frac(lattice, r)                        # 转回分数坐标用于周期插值
        V_line[idx] = sample_trilinear_periodic(V, f, grid) # LOCPOT 周期三线性插值
        rho_line[idx] = sample_trilinear_periodic(rho, f, grid)  # PARCHG 周期三线性插值
        if rho_chg_line is not None:                        # 如提供 CHGCAR
            rho_chg_line[idx] = sample_trilinear_periodic(rho_chg, f, grid)

    return s, V_line, rho_line, rho_chg_line                # 返回距离及对应剖面


def find_atoms_near_line(lattice, atom_fracs, atom_symbols, start_frac, direction_frac,
                         length, tol, periodic, highlight_tags):
    d_cart = frac_to_cart(lattice, direction_frac)
    d_norm = np.linalg.norm(d_cart)
    if d_norm < 1e-12:
        raise ValueError("Direction vector is zero in Cartesian space.")

    d_unit = d_cart / d_norm
    r0 = frac_to_cart(lattice, start_frac)

    if periodic:
        d_unit_frac = cart_to_frac(lattice, d_unit)
        end_frac = start_frac + d_unit_frac * length
        min_frac = np.minimum(start_frac, end_frac)
        max_frac = np.maximum(start_frac, end_frac)
        shift_min = np.floor(min_frac).astype(int) - 1
        shift_max = np.ceil(max_frac).astype(int) + 1
        shift_ranges = [
            range(shift_min[0], shift_max[0] + 1),
            range(shift_min[1], shift_max[1] + 1),
            range(shift_min[2], shift_max[2] + 1),
        ]
    else:
        shift_ranges = [range(0, 1), range(0, 1), range(0, 1)]
    markers = []
    for idx, (f0, sym) in enumerate(zip(atom_fracs, atom_symbols), start=1):
        best = None
        for sx in shift_ranges[0]:
            for sy in shift_ranges[1]:
                for sz in shift_ranges[2]:
                    f = f0 + np.array([sx, sy, sz], dtype=float)
                    r = frac_to_cart(lattice, f)
                    v = r - r0
                    s_proj = float(np.dot(v, d_unit))
                    s_clamped = min(max(s_proj, 0.0), length)
                    perp = v - s_clamped * d_unit
                    dist = np.linalg.norm(perp)
                    if dist <= tol + 1e-12:
                        if best is None or dist < best["dist"]:
                            best = {"s": s_clamped, "dist": dist}
        if best is not None:
            tag = highlight_tags.get(idx, "")
            markers.append({
                "s": best["s"],
                "element": sym,
                "index": idx,
                "tag": tag,
            })

    markers.sort(key=lambda m: m["s"])
    return markers


def build_element_colors(elements):
    unique = sorted(set(elements))
    cmap = plt.get_cmap("tab20" if len(unique) > 10 else "tab10")
    colors = {}
    for i, elem in enumerate(unique):
        colors[elem] = cmap(i % cmap.N)
    return colors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poscar", default="POSCAR", help="POSCAR (used for lattice reference)")
    ap.add_argument("--locpot", default="LOCPOT", help="LOCPOT file")
    ap.add_argument("--parchg", default="PARCHG", help="PARCHG file")
    ap.add_argument("--chgcar", default="CHGCAR", help="CHGCAR file (default: CHGCAR, auto-skip if missing)")
    ap.add_argument("--direction", nargs=3, type=float, required=True, metavar=("u", "v", "w"),
                    help="Fractional direction [u, v, w]")
    ap.add_argument("--start", nargs=3, type=float, required=True, metavar=("x0", "y0", "z0"),
                    help="Starting point [x0, y0, z0] in fractional coordinates")
    ap.add_argument("--length", type=float, default=None, help="Length (Å) along +direction (default: to boundary)")
    ap.add_argument("--step", type=float, default=None, help="Sampling step (Å) (default: VASP grid spacing)")
    ap.add_argument("--npts", type=int, default=None, help="Number of sampling points (override step)")
    ap.add_argument("--atom_tol", type=float, default=None, help="Atom-to-line tolerance (Å) (default: grid spacing)")
    ap.add_argument("--no_atom_periodic", action="store_true", help="Disable periodic images for atom marking")
    ap.add_argument("--out_prefix", default="line_profile", help="Output prefix")
    args = ap.parse_args()

    # Read lattice from POSCAR (authoritative geometry)
    with open(args.poscar, "r") as f:
        poscar_lines = f.readlines()
    lattice_poscar, *_ = read_poscar_header(poscar_lines)
    atom_fracs, atom_symbols = get_frac_positions_from_poscar(poscar_lines, lattice_poscar)

    # Read volumetric data
    lattice_L, grid_L, V = read_vasp_volumetric(args.locpot)
    lattice_C, grid_C, rho = read_vasp_volumetric(args.parchg)  # rho 指的是从 PARCHG 读出的三维网格数据
    rho_chg = None
    lattice_G = grid_G = None
    if args.chgcar:
        if os.path.exists(args.chgcar):
            lattice_G, grid_G, rho_chg = read_vasp_volumetric(args.chgcar)
        else:
            print(f"CHGCAR not found at {args.chgcar}, skip plotting CHGCAR.")

    # Basic consistency checks (grid must match for meaningful pairing)
    if grid_L != grid_C:
        raise RuntimeError(f"Grid mismatch: LOCPOT {grid_L} vs PARCHG {grid_C}")
    if rho_chg is not None and grid_L != grid_G:
        raise RuntimeError(f"Grid mismatch: LOCPOT {grid_L} vs CHGCAR {grid_G}")
    grid = grid_L

    # Use POSCAR lattice for coordinate conversions (recommended)
    lattice = lattice_poscar
    cell_volume = float(abs(np.linalg.det(lattice)))
    grid_volume = float(grid[0] * grid[1] * grid[2])

    # Convert PARCHG to density by dividing by cell volume (用户要求：rho = data / Vcell).
    rho = rho / cell_volume
    rho_unit = "Å^-3"
    rho_integral = float(rho.mean() * cell_volume)
    print(f"PARCHG converted to density ({rho_unit}); integral ≈ {rho_integral:.10g} e")

    # CHGCAR density（可选）
    rho_chg_line = None
    if rho_chg is not None:
        rho_chg = rho_chg / cell_volume
        rho_chg_integral = float(rho_chg.mean() * cell_volume)
        print(f"CHGCAR converted to density ({rho_unit}); integral ≈ {rho_chg_integral:.10g} e")

    direction = args.direction
    start_point = np.array(args.start, dtype=float)
    highlight_tags = {}

    if np.any(start_point < -1e-8) or np.any(start_point > 1 + 1e-8):
        print("Warning: start point outside [0,1], will be wrapped into the unit cell.")
        start_point = start_point - np.floor(start_point)

    s, V_line, rho_line, rho_chg_line = extract_profile_along_direction(
        lattice=lattice, grid=grid, V=V, rho=rho,
        direction_frac=direction, start_frac=start_point,
        length=args.length, step=args.step, npts=args.npts,
        rho_chg=rho_chg,
    )

    # Report potential range for quick inspection
    print(f"LOCPOT along line: V_min = {V_line.min():.6f} eV, V_max = {V_line.max():.6f} eV")

    # Geometry visualization for sanity check
    geom_png = f"{args.out_prefix}_geom.png"
    try:
        plot_line_geometry(lattice, start_point, direction, float(s[-1]), geom_png)
        print("Saved geometry figure:", geom_png)
    except Exception as e:
        print("Geometry plot failed:", e)

    # Save DAT
    out_dat = f"{args.out_prefix}.dat"
    if rho_chg_line is not None:
        header = f"s(Angstrom)\tV_locpot(eV)\trho_parchg({rho_unit})\trho_chgcar({rho_unit})"
        data_out = np.column_stack([s, V_line, rho_line, rho_chg_line])
    else:
        header = f"s(Angstrom)\tV_locpot(eV)\trho_parchg({rho_unit})"
        data_out = np.column_stack([s, V_line, rho_line])
    np.savetxt(out_dat, data_out, header=header, delimiter="\t", fmt="%.10e")

    # Plot (Fig4-like: potential + density vs distance)
    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    color_v = "tab:blue"
    color_rho = "tab:orange"
    ax1.plot(s, V_line, label="V (LOCPOT)", color=color_v)
    ax1.set_xlabel("s (Å) along direction")
    ax1.set_ylabel("V (eV)", color=color_v)
    ax1.tick_params(axis="y", labelcolor=color_v)
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(s, rho_line, label="rho (PARCHG)", color=color_rho)
    if rho_chg_line is not None:
        ax2.plot(s, rho_chg_line, label="rho (CHGCAR)", color="tab:green")
    ax2.set_ylabel(f"Density ({rho_unit})", color=color_rho)
    ax2.tick_params(axis="y", labelcolor=color_rho)
    
    atom_tol = args.atom_tol
    if atom_tol is None:
        atom_tol = compute_grid_step(lattice, grid)
    markers = find_atoms_near_line(
        lattice=lattice,
        atom_fracs=atom_fracs,
        atom_symbols=atom_symbols,
        start_frac=start_point,
        direction_frac=direction,
        length=float(s[-1]),
        tol=atom_tol,
        periodic=not args.no_atom_periodic,
        highlight_tags=highlight_tags,
    )

    element_colors = None
    if markers:
        element_colors = build_element_colors([m["element"] for m in markers])
        label_dx = compute_grid_step(lattice, grid)
        used_x = []
        for m in markers:
            x = m["s"]
            color = element_colors[m["element"]]
            lw = 1.5 if m["tag"] else 1.0
            ax1.axvline(x, color=color, linestyle="--", linewidth=lw, alpha=0.7)
            level = sum(1 for ux in used_x if abs(ux - x) < label_dx)
            y = 1.02 + 0.05 * min(level, 6)
            label = f"{m['element']}{m['index']}{m['tag']}"
            ax1.text(
                x, y, label,
                transform=ax1.get_xaxis_transform(),
                color=color,
                ha="center",
                va="bottom",
                fontsize=8,
                fontweight="bold" if m["tag"] else "normal",
            )
            used_x.append(x)

    # Legend merge
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    element_handles = []
    element_labels = []
    if element_colors is not None:
        for elem, color in element_colors.items():
            element_handles.append(Line2D([0], [0], color=color, linestyle="--"))
            element_labels.append(elem)

    ax1.legend(h1 + h2 + element_handles, l1 + l2 + element_labels, loc="best")

    out_png = f"{args.out_prefix}.png"
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print("Done.")
    print("Saved:", out_dat)
    print("Saved:", out_png)


if __name__ == "__main__":
    main()
