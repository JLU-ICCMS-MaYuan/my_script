#!/usr/bin/env python3
import os
import re
import glob
import argparse
import matplotlib.pyplot as plt

def natural_sort_key(s):
    """Sort strings containing numbers in a natural way (e.g., 2 before 10)."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def get_all_gibbs(filename):
    values = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if "Gibbs Free energy =" in line:
                    match = re.search(r"Gibbs Free energy\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
                    if match:
                        values.append(float(match.group(1)))
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return values

def main():
    parser = argparse.ArgumentParser(description="Extract Gibbs energy from SSCHA output files.")
    parser.add_argument("-n", "--natoms", type=int, default=1, help="Total number of atoms (default: 1)")
    args = parser.parse_args()

    # Get all sscha_*.out files and sort them naturally
    files = glob.glob("sscha_*.out")
    files.sort(key=natural_sort_key)
    
    if not files:
        print("未找到 sscha_*.out 文件。")
        return

    all_data = []
    for f in files:
        vals = get_all_gibbs(f)
        if vals:
            for i, v in enumerate(vals):
                all_data.append({
                    "file": f,
                    "index": i + 1,
                    "value": v
                })
        else:
            print(f"警告: 在 {f} 中未找到 'Gibbs Free energy'")

    if not all_data:
        print("在任何文件中都未找到 Gibbs 能量值。")
        return

    # 以最后一个提取到的值为基准
    reference = all_data[-1]
    ref_val = reference["value"]
    
    print(f"\n参数: 总原子数 = {args.natoms}")
    print(f"以 {reference['file']} 的第 {reference['index']} 个 Gibbs 值为基准 (Gibbs = {ref_val:.10e} eV)\n")
    print(f"{ '文件名':<20} | { '序号':<6} | { 'Gibbs (eV)':<20} | { '差值 (meV/atom)':<20}")
    print("-" * 85)
    
    diffs_meV_atom = []
    labels = []
    
    for item in all_data:
        # 差值 (meV/atom) = (val - ref_val) * 1000 / natoms
        diff = (item["value"] - ref_val) * 1000.0 / args.natoms
        diffs_meV_atom.append(diff)
        
        # 为绘图创建标签
        label = f"{item['file']}({item['index']})"
        labels.append(label)
        
        print(f"{item['file']:<20} | {item['index']:<6} | {item['value']:<20.6f} | {diff:<20.6f}")

    # 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(diffs_meV_atom)), diffs_meV_atom, marker='o', linestyle='-', color='b')
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
    plt.ylabel("Gibbs Energy Difference (meV/atom)")
    plt.xlabel("Output Step")
    plt.title(f"Gibbs Energy Convergence (Reference: {reference['file']} step {reference['index']})")
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    plt.savefig("Gibbs_energy.png")
    print(f"\n图像已保存为: Gibbs_energy.png")

if __name__ == "__main__":
    main()
