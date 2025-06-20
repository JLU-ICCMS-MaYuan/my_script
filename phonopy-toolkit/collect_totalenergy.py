#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from typing import Optional

# ---------- 解析函数 ---------- #
def get_num_atoms_vasp(outcar_path: str) -> Optional[int]:
    """从 VASP OUTCAR 里解析原子数（出现错误返回 None）"""
    try:
        res = subprocess.run(
            r'grep -n -s "position of ions in cartesian coordinates" ' + outcar_path,
            shell=True, capture_output=True, text=True,
        )
        if not res.stdout:
            return None
        begin_id = int(res.stdout.split(":")[0])
        natoms = 0
        with open(outcar_path) as f:
            for i, line in enumerate(f, 1):
                if i > begin_id:
                    coords = re.findall(r"[-+]?\d+\.\d+", line)
                    if len(coords) == 3:
                        natoms += 1
                    else:
                        break
        return natoms or None
    except Exception:
        return None


def get_num_atoms_qe(scf_path: str) -> Optional[int]:
    """从 QE scf.out 里解析原子数（出现错误返回 None）"""
    try:
        with open(scf_path) as f:
            for line in f:
                if "number of atoms/cell" in line:
                    # 例如：  number of atoms/cell =  20
                    m = re.search(r"=\s*([0-9]+)", line)
                    if m:
                        return int(m.group(1))
        return None
    except Exception:
        return None


def get_energy_vasp(outcar_path: str) -> Optional[float]:
    """最后一次 free‑energy TOTEN (eV)"""
    try:
        res = subprocess.run(
            r'grep -s "free  energy   TOTEN" ' + outcar_path + " | tail -1 | awk '{print $5}'",
            shell=True, capture_output=True, text=True,
        )
        return float(res.stdout.strip()) if res.stdout.strip() else None
    except Exception:
        return None


def get_energy_qe(scf_path: str) -> Optional[float]:
    """最后一次 “!    total energy” 行的能量 (Ry -> eV)"""
    try:
        # grep, 取最后一行
        res = subprocess.run(
            r'grep -s "!    total energy" ' + scf_path + " | tail -1 | awk '{print $5}'",
            shell=True, capture_output=True, text=True,
        )
        en_ry = res.stdout.strip()
        return float(en_ry) * 13.6056980659 if en_ry else None  # Ry→eV
    except Exception:
        return None


# ---------- 振幅生成 ---------- #
def generate_amplitudes(start: float, stop: float, step: float) -> list[float]:
    a = []
    cur = start
    # 解决浮点误差，保证包含 stop
    while (step > 0 and cur <= stop + 1e-12) or (step < 0 and cur >= stop - 1e-12):
        a.append(round(cur, 10))  # 避免累积误差
        cur += step
    return a


# ---------- 主程序 ---------- #
def main():
    parser = argparse.ArgumentParser(
        description="Collect energies from MPOSCAR_* directories."
    )
    parser.add_argument("-amp", "--amplitudes",
                        nargs=3, type=float, default=[-3.0, 3.0, 0.1],
                        help="start stop step (default: -3.0 3.0 0.1)")
    parser.add_argument("--code", choices=["vasp", "qe", "auto"],
                        default="auto",
                        help="Which code to parse: vasp | qe | auto (default auto)")
    args = parser.parse_args()

    amps = generate_amplitudes(*args.amplitudes)
    print(f"{'Amplitude':>10}  {'Natoms':>6}  {'Energy (eV/atom)':>18}  {'Source':>8}")

    for amp in amps:
        d = f"MPOSCAR_{amp:0.2f}"
        outcar  = os.path.join(d, "OUTCAR")
        scfout  = os.path.join(d, "scf.out")

        energy, natoms, src = 1e14, 1, "None"   # 默认失败

        if args.code in ("vasp", "auto") and os.path.isfile(outcar):
            e = get_energy_vasp(outcar)
            n = get_num_atoms_vasp(outcar)
            if e is not None:
                energy, src = e, "VASP"
            if n is not None:
                natoms = n

        if args.code in ("qe", "auto") and src == "None" and os.path.isfile(scfout):
            e = get_energy_qe(scfout)
            n = get_num_atoms_qe(scfout)
            if e is not None:
                energy, src = e, "QE"
            if n is not None:
                natoms = n

        print(f"{amp:10.2f}  {natoms:6d}  {energy/natoms:18.8f}  {src:>8}")


if __name__ == "__main__":
    main()
