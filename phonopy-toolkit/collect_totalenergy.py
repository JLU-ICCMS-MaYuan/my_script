#!/usr/bin/env python3
import argparse
import subprocess
import os
import re

def get_num_atoms(outcar_path):
    try:
        result = subprocess.run(f'grep -n -s "position of ions in cartesian coordinates" {outcar_path}', 
                                shell=True, capture_output=True, text=True)
        begin_id = int(result.stdout.split(":")[0])

        N = 0
        with open(outcar_path, 'r') as f:
            for i, line in enumerate(f, start=1):
                if i > begin_id:
                    coords = re.findall(r"[-+]?\d+\.\d+", line.strip())
                    if len(coords) == 3:
                        N += 1
                    else:
                        break
        return N
    except Exception:
        return 1


def get_free_energy(outcar_path):
    """Extract the free energy (TOTEN) from an OUTCAR file."""
    try:
        result = subprocess.run(
            f'grep -s "free  energy   TOTEN" {outcar_path} | tail -n 1 | awk \'{{print $5}}\'',
            shell=True, capture_output=True, text=True
        )
        dE = result.stdout.strip()
        return float(dE) if dE else 1e14  # Return a large number if no energy is found
    except Exception:
        return 1e14  # Return a large number if an error occurs


def generate_amplitudes(start: float, stop: float, step: float) -> list[float]:
    """Generate a list of amplitudes based on start, stop, and step."""
    amplitudes = []
    current = start
    while current <= stop:
        amplitudes.append(current)
        current += step
    return amplitudes


def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Generate force constants and modulation for a given POSCAR.")
    parser.add_argument(
        "-amp",
        "--amplitudes",
        nargs=3,
        type=float,
        default=[-3.0, 3.0, 0.1],
        help="Three values: start, stop, and step for amplitudes (default: -3.0 3.0 0.1).",
    )
    args = parser.parse_args()

    # Generate amplitudes based on start, stop, and step
    amplitudes = generate_amplitudes(args.amplitudes[0], args.amplitudes[1], args.amplitudes[2])

    # Scan the current directory for OUTCAR files and extract free energy
    outcar_files = [f for f in os.listdir() if f.startswith("OUTCAR")]
    print("{:>30}, {:>9}, {:>6}, {:>15}".format("OUTCAR_file", "Amplitude", "natoms", "total Energy"))
    for idx, amp in enumerate(amplitudes):
        # outcar_file = os.path.join(f"{idx+1}.modu_delta{amp:.2f}", "OUTCAR")
        outcar_file = os.path.join(f"modu_delta{amp:.2f}", "OUTCAR")
        energy = get_free_energy(outcar_file)
        natoms = get_num_atoms(outcar_file)
        print("{:>30}, {:>9.2f}, {:>6}, {:>15.8f}".format(outcar_file, amp, natoms, energy/natoms))


if __name__ == "__main__":
    main()
