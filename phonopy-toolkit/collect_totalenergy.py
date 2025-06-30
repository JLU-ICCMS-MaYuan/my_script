#!/usr/bin/env python3
import argparse
import os
import re
import glob
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
    try:
        # grep, 取最后一行
        res = subprocess.run(
            r'grep -s "!    total energy" ' + scf_path + " | tail -1 | awk '{print $5}'",
            shell=True, capture_output=True, text=True,
        )
        en_ry = res.stdout.strip()
        return float(en_ry) if en_ry else None  # Ry→eV
    except Exception:
        return None


# ---------- 振幅生成 ---------- #
def generate_amplitudes(start: float, stop: float, step: float) -> list[float]:
    if start == 0. and stop == 0. and step == 0.:
        a = []
        # 匹配所有 MPOSCAR_* 开头的目录
        for path in glob.glob("MPOSCAR_*"):
            if os.path.isdir(path):
                try:
                    val = float(path.split("_")[1])
                    a.append(val)
                except (IndexError, ValueError):
                    continue
        return sorted(a)
    
    elif start != 0. and stop != 0. and step != 0.:
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
    parser.add_argument("-a", "--amplitudes",
                        nargs=3, type=float, default=[0., 0., 0.],
                        help="start stop step (default: -3.0 3.0 0.1)")
    parser.add_argument("-c", "--code", choices=["vasp", "qe"],
                        default="vasp",
                        help="Which code to parse: vasp | qe ")
    args = parser.parse_args()

    amps = generate_amplitudes(*args.amplitudes)

    if args.code == "qe":
        print(f"{'Amplitude':>10}    {'Energy (Ry/atom)':>18}")
        for amp in amps:
            d = f"MPOSCAR_{amp:0.2f}"
            scfout  = os.path.join(d, "scf.out")
            energy, natoms = 1e14, 1
            if os.path.isfile(scfout):
                energy = get_energy_qe(scfout)
                natoms = get_num_atoms_qe(scfout)
            print(f"{amp:10.2f} {energy/natoms:18.8f}")
    elif args.code == "vasp":
        print(f"{'Amplitude':>10}    {'Energy (eV/atom)':>18}")
        for amp in amps:
            d = f"MPOSCAR_{amp:0.2f}"
            outcar  = os.path.join(d, "scf/OUTCAR")
            energy, natoms = 1e14, 1
            if os.path.isfile(outcar):
                energy = get_energy_vasp(outcar)
                natoms = get_num_atoms_vasp(outcar)
            print(f"{amp:10.2f} {energy/natoms:18.8f}")
    else:
        print("You have to specify qe or vasp")


if __name__ == "__main__":
    main()
