#!/usr/bin/env python3
"""
XSF 压强检查脚本：直接从 XSF 文件中提取压强信息，分析偏离目标压强的结构

使用方法：在 XSF 目录下运行此脚本
    cd XSF
    python3 check_xsf_pressure.py -p 100

压强计算原理（基于 vasp2xsf.py 转换逻辑）：
    - VIRIAL = -stress × volume
    - VIRIAL 单位：eV（电子伏特）
    - 压强计算：stress = -VIRIAL / volume（单位：eV/Å³）
    - 转换为 GPa：1 eV/Å³ = 160.21766208 GPa
"""

import os
import re
import argparse
import glob
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import shutil
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


def extract_pressure_from_xsf(xsf_path):
    """
    从 XSF 文件中提取压强信息

    主要方法：从 VIRIAL 矩阵计算压强
    VIRIAL 格式：
        VIRIAL
                 89.18844723  1.43016552  1.65702783  1.43016552  81.15942903  2.68153854  1.65702783  2.68153854  96.56524834

    这是一个 3x3 矩阵写成一行，对角元（第 0, 4, 8 个数）的平均值即为压强

    单位关系（基于 vasp2xsf.py 转换脚本）：
        VIRIAL = -stress × volume
        其中 stress 单位是 eV/Å³，volume 单位是 Å³
        所以 VIRIAL 单位是 eV（电子伏特），不是 eV/Å³

    计算压强：
        1. 读取 VIRIAL 对角元（单位：eV）
        2. 读取晶胞参数 PRIMVEC，计算体积（单位：Å³）
        3. 计算 stress = -VIRIAL / volume（单位：eV/Å³）
        4. 计算 pressure = (stress_xx + stress_yy + stress_zz) / 3（单位：eV/Å³）
        5. 转换为 GPa：1 eV/Å³ = 160.21766208 GPa

    备用方法：如果没有 VIRIAL，尝试从注释中读取

    返回: (total_pressure_gpa, source_info) 如果未找到则返回 (None, None)
    """
    total_pressure_gpa = None
    source_info = None

    try:
        with open(xsf_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

            # 方法1: 从 VIRIAL 矩阵读取（优先）
            virial_data = None
            primvec_data = None

            # 先找 VIRIAL 和 PRIMVEC
            for i, line in enumerate(lines):
                if line.strip() == 'VIRIAL':
                    # 读取下一行的 9 个数字
                    if i + 1 < len(lines):
                        next_line = lines[i + 1].strip()
                        numbers = re.findall(r'[-\d.]+', next_line)
                        if len(numbers) >= 9:
                            virial_data = [float(n) for n in numbers[:9]]

                if line.strip() == 'PRIMVEC':
                    # 读取接下来的3行晶胞矢量
                    if i + 3 < len(lines):
                        primvec_lines = []
                        for j in range(1, 4):
                            vec_line = lines[i + j].strip()
                            vec_numbers = re.findall(r'[-\d.]+', vec_line)
                            if len(vec_numbers) >= 3:
                                primvec_lines.append([float(n) for n in vec_numbers[:3]])
                        if len(primvec_lines) == 3:
                            primvec_data = np.array(primvec_lines)

            # 如果找到了 VIRIAL 和 PRIMVEC，计算压强
            if virial_data is not None and primvec_data is not None:
                # 取 VIRIAL 对角元（index 0, 4, 8）
                virial_xx = virial_data[0]  # eV
                virial_yy = virial_data[4]  # eV
                virial_zz = virial_data[8]  # eV

                # 计算晶胞体积（Å³）
                # volume = a · (b × c)
                volume = np.abs(np.dot(primvec_data[0], np.cross(primvec_data[1], primvec_data[2])))

                # 计算应力：stress = -VIRIAL / volume（单位：eV/Å³）
                stress_xx = -virial_xx / volume
                stress_yy = -virial_yy / volume
                stress_zz = -virial_zz / volume

                # 静水压强：对角元平均值（单位：eV/Å³）
                pressure_ev_per_a3 = (stress_xx + stress_yy + stress_zz) / 3.0

                # 转换为 GPa：1 eV/Å³ = 160.21766208 GPa
                total_pressure_gpa = pressure_ev_per_a3 * 160.21766208

                source_info = f"VIRIAL/volume ({pressure_ev_per_a3:.6f} eV/Å³, V={volume:.2f} Å³)"

                return total_pressure_gpa, source_info

            # 方法2: 从注释中读取（备用）
            for line in lines:
                # 格式1: # Pressure: 123.45 GPa
                match = re.search(r'#\s*Pressure\s*[=:]\s*([-\d.]+)\s*GPa', line, re.IGNORECASE)
                if match:
                    total_pressure_gpa = float(match.group(1))
                    source_info = "Pressure comment"
                    return total_pressure_gpa, source_info

                # 格式2: # Total Pressure = 123.45 GPa
                match = re.search(r'#\s*Total\s+Pressure\s*[=:]\s*([-\d.]+)\s*GPa', line, re.IGNORECASE)
                if match:
                    total_pressure_gpa = float(match.group(1))
                    source_info = "Total Pressure comment"
                    return total_pressure_gpa, source_info

                # 格式3: # external pressure = 1234.5 kB  Pullay stress = 10000.0 kB
                match = re.search(r'#\s*external pressure\s*=\s*([-\d.]+)\s*kB.*?Pullay stress\s*=\s*([-\d.]+)\s*kB', line, re.IGNORECASE)
                if match:
                    ext_press_kb = float(match.group(1))
                    pullay_kb = float(match.group(2))
                    total_pressure_gpa = (ext_press_kb + pullay_kb) / 10.0
                    source_info = "external + Pullay comment"
                    return total_pressure_gpa, source_info

    except Exception as e:
        print(f"[WARNING] 无法读取 {xsf_path}: {e}")

    return total_pressure_gpa, source_info


def extract_it_number(xsf_path):
    """从 XSF 文件路径中提取 IT 编号，例如 IT9/IT9_SCF_structure.xsf -> 9"""
    match = re.search(r'IT(\d+)', str(xsf_path))
    return int(match.group(1)) if match else -1


def analyze_xsf_pressures(target_pressure_gpa, threshold_percent=20.0, extreme_error_percent=None, move_extreme_count=None):
    """
    分析所有 XSF 文件中的压强信息

    参数:
        target_pressure_gpa: 目标压强 (GPa)
        threshold_percent: 允许的合理阈值 (百分比)
        extreme_error_percent: 极端不合理的最小阈值 (百分比)
        move_extreme_count: 移动极端不合理结构的前 n 个到垃圾箱（None 表示不移动）
    """
    # 检查当前目录
    cwd = Path.cwd()
    if cwd.name != "XSF":
        print(f"[ERROR] 当前目录不是 XSF 目录")
        print(f"[ERROR] 当前目录: {cwd}")
        print(f"[ERROR] 请先切换到 XSF 目录: cd XSF")
        return

    # 计算合理阈值范围（基于百分比）
    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    print(f"[INFO] 当前工作目录: {cwd}")
    print(f"[INFO] 目标压强: {target_pressure_gpa} GPa")
    print(f"[INFO] VIRIAL 单位: eV (电子伏特)")
    print(f"[INFO] 压强计算: stress = -VIRIAL / volume (单位: eV/Å³ → GPa)")
    print(f"[INFO] 合理阈值: ±{threshold_percent}%")
    print(f"[INFO] 合理范围: {threshold_lower:.2f} ~ {threshold_upper:.2f} GPa")

    # 计算极端不合理范围（基于百分比）
    if extreme_error_percent is not None:
        extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
        extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)
        print(f"[INFO] 极端不合理阈值: ±{extreme_error_percent}%")
        print(f"[INFO] 极端不合理范围: < {extreme_lower:.2f} GPa 或 > {extreme_upper:.2f} GPa")
        print(f"[INFO] 超出此范围的结构将不在主图中显示，单独处理")
    print()

    # 查找所有 XSF 文件（IT*/IT*_SCF_*.xsf 格式）
    xsf_pattern = "IT*/IT*_SCF_*.xsf"
    xsf_files = glob.glob(xsf_pattern, recursive=False)

    if not xsf_files:
        print(f"[ERROR] 未找到任何匹配的 XSF 文件: {xsf_pattern}")
        print(f"[ERROR] 请确认：")
        print(f"[ERROR]   1. 当前在 XSF 目录下")
        print(f"[ERROR]   2. 存在 IT*/IT*_SCF_*.xsf 文件")
        return

    print(f"[INFO] 找到 {len(xsf_files)} 个 XSF 文件")
    print()

    # 存储所有压强数据
    pressure_data = []  # 正常数据（在极大误差范围内，用于主图）
    outliers = []  # 偏离目标的结构（超过threshold但在extreme范围内）
    extreme_errors = []  # 极大误差的结构（超过extreme范围）
    no_pressure_files = []  # 没有压强信息的文件

    # 提取压强信息（并行处理）
    print("[INFO] 正在提取 XSF 文件中的压强信息...")

    # 确定使用的进程数（最多使用 CPU 核心数，但至少保留1个核心）
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"[INFO] 使用 {num_workers} 个进程并行处理")

    # 使用进程池并行处理
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # 提交所有任务
        future_to_path = {executor.submit(extract_pressure_from_xsf, path): path for path in xsf_files}

        # 使用 tqdm 显示进度
        for future in tqdm(as_completed(future_to_path), total=len(xsf_files), desc="处理 XSF", unit="文件"):
            xsf_path = future_to_path[future]
            try:
                total_pressure_gpa, source_info = future.result()
            except Exception as e:
                tqdm.write(f"[WARNING] 处理 {xsf_path} 时出错: {e}")
                continue

            if total_pressure_gpa is not None:
                data_entry = {
                    'xsf_path': xsf_path,
                    'total_pressure_gpa': total_pressure_gpa,
                    'source_info': source_info,
                    'it_number': extract_it_number(xsf_path)
                }

                # 检查是否为极大误差
                is_extreme = False
                if extreme_error_percent is not None:
                    extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
                    extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)
                    if total_pressure_gpa < extreme_lower or total_pressure_gpa > extreme_upper:
                        is_extreme = True
                        extreme_errors.append({
                            'xsf_path': xsf_path,
                            'total_pressure_gpa': total_pressure_gpa,
                            'deviation_gpa': abs(total_pressure_gpa - target_pressure_gpa),
                            'it_number': data_entry['it_number']
                        })

                # 如果不是极大误差，加入主数据
                if not is_extreme:
                    pressure_data.append(data_entry)

                    # 检查是否偏离合理范围
                    if total_pressure_gpa < threshold_lower or total_pressure_gpa > threshold_upper:
                        outliers.append({
                            'xsf_path': xsf_path,
                            'total_pressure_gpa': total_pressure_gpa,
                            'deviation_gpa': abs(total_pressure_gpa - target_pressure_gpa),
                            'it_number': data_entry['it_number']
                        })
            else:
                no_pressure_files.append(xsf_path)

    if not pressure_data and not extreme_errors:
        print("[ERROR] 未能从任何 XSF 文件中提取到压强信息")
        if no_pressure_files:
            print(f"[ERROR] 有 {len(no_pressure_files)} 个文件不包含压强信息")
            print("[INFO] 请确认 XSF 文件中是否包含压强信息（通常在注释中）")
        return

    print(f"[INFO] 成功提取 {len(pressure_data) + len(extreme_errors)} 个结构的压强信息")
    if no_pressure_files:
        print(f"[WARNING] 有 {len(no_pressure_files)} 个文件不包含压强信息")
    if extreme_error_percent is not None:
        print(f"[INFO]   - 主图数据: {len(pressure_data)} 个结构")
        print(f"[INFO]   - 极端不合理: {len(extreme_errors)} 个结构 (超过 ±{extreme_error_percent}%)")
    print(f"[INFO] 发现 {len(outliers)} 个偏离合理范围（±{threshold_percent}%）的结构（主图中）")
    print()

    # 统计信息（主图数据）
    if pressure_data:
        total_pressures = [d['total_pressure_gpa'] for d in pressure_data]
        print("=" * 70)
        print("压强统计信息 (主图数据，GPa)")
        print("=" * 70)
        print(f"  最小值: {min(total_pressures):.2f} GPa")
        print(f"  最大值: {max(total_pressures):.2f} GPa")
        print(f"  平均值: {np.mean(total_pressures):.2f} GPa")
        print(f"  中位数: {np.median(total_pressures):.2f} GPa")
        print(f"  标准差: {np.std(total_pressures):.2f} GPa")
        print("=" * 70)
        print()

        # 绘制压强分布图
        plot_pressure_distribution(pressure_data, target_pressure_gpa, threshold_percent)

    # 处理偏离的结构
    if outliers:
        save_outliers_to_file(outliers, target_pressure_gpa, threshold_percent, threshold_lower, threshold_upper)

    # 处理极大误差的结构
    if extreme_errors:
        save_extreme_errors_to_file(extreme_errors, target_pressure_gpa, extreme_error_percent)

        # 如果指定了移动参数，移动前 n 个极端不合理的结构
        if move_extreme_count is not None and move_extreme_count > 0:
            structures_to_move = extreme_errors[:move_extreme_count]
            move_extreme_errors_to_trash(structures_to_move, move_extreme_count, len(extreme_errors))

        # 绘制极大误差结构的分布图
        plot_extreme_error_distribution(extreme_errors, target_pressure_gpa, extreme_error_percent)


def save_outliers_to_file(outliers, target_pressure_gpa, threshold_percent, threshold_lower, threshold_upper):
    """保存偏离结构列表到文件"""
    print("=" * 70)
    print(f"偏离合理范围（±{threshold_percent}%）的结构列表")
    print("=" * 70)

    # 按偏差从大到小排序
    outliers.sort(key=lambda x: x['deviation_gpa'], reverse=True)

    # 检测文件是否存在
    missing_count = 0
    for outlier in outliers:
        xsf_path = outlier['xsf_path']
        outlier['file_exists'] = Path(xsf_path).exists()
        if not outlier['file_exists']:
            missing_count += 1

    # 输出到文件
    with open('outlier_xsf_structures.txt', 'w', encoding='utf-8') as f:
        f.write(f"# 偏离合理范围的 XSF 结构 (目标压强 {target_pressure_gpa} GPa, 合理阈值 ±{threshold_percent}%)\n")
        f.write(f"# 共 {len(outliers)} 个结构\n")
        f.write(f"# 合理范围: {threshold_lower:.2f} ~ {threshold_upper:.2f} GPa\n")
        if missing_count > 0:
            f.write(f"# ⚠️ 警告: 有 {missing_count} 个 XSF 文件不存在（标记为 [文件不存在]）\n")
        f.write("#\n")
        f.write("# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa) | IT编号 | [状态]\n")
        f.write("#" + "=" * 68 + "\n\n")

        for i, outlier in enumerate(outliers, 1):
            xsf_path = outlier['xsf_path']
            pressure = outlier['total_pressure_gpa']
            deviation = outlier['deviation_gpa']
            it_num = outlier['it_number']
            file_exists = outlier.get('file_exists', False)

            status = "" if file_exists else "  [文件不存在]"
            line = f"{i:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa  |  IT{it_num}{status}\n"
            print(f"  {line.strip()}")
            f.write(line)

    print()
    print(f"[INFO] 偏离结构列表已保存到: outlier_xsf_structures.txt")
    if missing_count > 0:
        print(f"[WARNING] ⚠️ 有 {missing_count} 个 XSF 文件不存在，已在文件中标记")
    print()


def save_extreme_errors_to_file(extreme_errors, target_pressure_gpa, extreme_error_percent):
    """保存极大误差结构列表到文件"""
    print("=" * 70)
    print(f"极大误差结构列表 (超过目标压强 ±{extreme_error_percent}%)")
    print("=" * 70)

    # 统计极大误差结构的压强信息
    extreme_pressures = [e['total_pressure_gpa'] for e in extreme_errors]
    print()
    print("极大误差结构压强统计 (GPa)")
    print("-" * 70)
    print(f"  最小值: {min(extreme_pressures):.2f} GPa")
    print(f"  最大值: {max(extreme_pressures):.2f} GPa")
    print(f"  平均值: {np.mean(extreme_pressures):.2f} GPa")
    print(f"  中位数: {np.median(extreme_pressures):.2f} GPa")
    print(f"  标准差: {np.std(extreme_pressures):.2f} GPa")
    print("-" * 70)
    print()

    # 按偏差从大到小排序
    extreme_errors.sort(key=lambda x: x['deviation_gpa'], reverse=True)

    # 检测文件是否存在
    missing_count = 0
    for error in extreme_errors:
        xsf_path = error['xsf_path']
        error['file_exists'] = Path(xsf_path).exists()
        if not error['file_exists']:
            missing_count += 1

    # 输出到文件
    with open('extreme_error_xsf_structures.txt', 'w', encoding='utf-8') as f:
        extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
        extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)

        f.write(f"# 极大误差 XSF 结构列表 (超过目标压强 {target_pressure_gpa} GPa 的 ±{extreme_error_percent}%)\n")
        f.write(f"# 共 {len(extreme_errors)} 个结构\n")
        f.write(f"# 极大误差范围: < {extreme_lower:.2f} GPa 或 > {extreme_upper:.2f} GPa\n")
        if missing_count > 0:
            f.write(f"# ⚠️ 警告: 有 {missing_count} 个 XSF 文件不存在（标记为 [文件不存在]）\n")
        f.write("#\n")
        f.write("# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa) | IT编号 | [状态]\n")
        f.write("#" + "=" * 68 + "\n\n")

        for i, error in enumerate(extreme_errors, 1):
            xsf_path = error['xsf_path']
            pressure = error['total_pressure_gpa']
            deviation = error['deviation_gpa']
            it_num = error['it_number']
            file_exists = error.get('file_exists', False)

            status = "" if file_exists else "  [文件不存在]"
            line = f"{i:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa  |  IT{it_num}{status}\n"
            if i <= 10:  # 只打印前10个到终端
                print(f"  {line.strip()}")
            f.write(line)

        if len(extreme_errors) > 10:
            print(f"  ... (还有 {len(extreme_errors) - 10} 个，详见文件)")

    print()
    print(f"[INFO] 极大误差结构列表已保存到: extreme_error_xsf_structures.txt")
    if missing_count > 0:
        print(f"[WARNING] ⚠️ 有 {missing_count} 个 XSF 文件不存在，已在文件中标记")
    print()


def plot_pressure_distribution(pressure_data, target_pressure_gpa, threshold_percent):
    """绘制压强分布图"""
    total_pressures = [d['total_pressure_gpa'] for d in pressure_data]

    # 计算合理范围
    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # 子图1: 直方图
    n, bins, patches = ax1.hist(total_pressures, bins=50, alpha=0.7, color='steelblue', edgecolor='black')

    # 填充合理范围
    ax1.axvspan(threshold_lower, threshold_upper,
                alpha=0.15, color='green', label=f'Reasonable Range (±{threshold_percent}%)', zorder=1)

    # 标记目标压强
    ax1.axvline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=3)
    ax1.axvline(threshold_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Upper Limit: {threshold_upper:.2f} GPa', zorder=2)
    ax1.axvline(threshold_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Lower Limit: {threshold_lower:.2f} GPa', zorder=2)

    # 添加文本标注
    y_max = n.max() * 0.95
    ax1.text(target_pressure_gpa, y_max, f'{target_pressure_gpa} GPa',
             ha='center', va='top', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))

    ax1.set_xlabel('Total Pressure (GPa)', fontsize=12)
    ax1.set_ylabel('Number of Structures', fontsize=12)
    ax1.set_title('Pressure Distribution of XSF Structures (Histogram)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)

    # 子图2: 散点图（按序号排列）
    indices = list(range(1, len(total_pressures) + 1))
    colors = ['green' if threshold_lower <= p <= threshold_upper else 'red' for p in total_pressures]

    ax2.axhspan(threshold_lower, threshold_upper,
                alpha=0.15, color='green', label=f'Reasonable Range (±{threshold_percent}%)', zorder=1)

    ax2.scatter(indices, total_pressures, c=colors, alpha=0.6, s=20, zorder=3)

    ax2.axhline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=2)
    ax2.axhline(threshold_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Upper Limit: {threshold_upper:.2f} GPa', zorder=2)
    ax2.axhline(threshold_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Lower Limit: {threshold_lower:.2f} GPa', zorder=2)

    ax2.set_xlabel('Structure Index', fontsize=12)
    ax2.set_ylabel('Total Pressure (GPa)', fontsize=12)
    ax2.set_title('Pressure Distribution of XSF Structures (Scatter Plot)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    output_file = 'xsf_pressure_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"[INFO] 压强分布图已保存到: {output_file}")
    print()


def plot_extreme_error_distribution(extreme_errors, target_pressure_gpa, extreme_error_percent):
    """绘制极大误差结构的压强分布图（带 IT 分隔线）"""
    extreme_pressures = [e['total_pressure_gpa'] for e in extreme_errors]

    # 按 IT 编号排序
    extreme_errors_sorted_by_it = sorted(extreme_errors, key=lambda x: x['it_number'])

    # 找出 IT 边界
    it_boundaries = []
    current_it = None
    for idx, error in enumerate(extreme_errors_sorted_by_it):
        if current_it is not None and error['it_number'] != current_it:
            it_boundaries.append(idx)
        current_it = error['it_number']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # 计算极大误差范围
    extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
    extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)

    # 子图1: 直方图
    n, bins, patches = ax1.hist(extreme_pressures, bins=30, alpha=0.7, color='crimson', edgecolor='black')

    ax1.axvline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=3)
    ax1.axvline(extreme_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Upper: {extreme_upper:.2f} GPa', zorder=2)
    ax1.axvline(extreme_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Lower: {extreme_lower:.2f} GPa', zorder=2)

    ax1.axvspan(ax1.get_xlim()[0], extreme_lower, alpha=0.15, color='red', label='Extreme Error Zone', zorder=1)
    ax1.axvspan(extreme_upper, ax1.get_xlim()[1], alpha=0.15, color='red', zorder=1)

    y_max = n.max() * 0.95
    ax1.text(target_pressure_gpa, y_max, f'{target_pressure_gpa} GPa',
             ha='center', va='top', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))

    ax1.set_xlabel('Total Pressure (GPa)', fontsize=12)
    ax1.set_ylabel('Number of Structures', fontsize=12)
    ax1.set_title(f'Extreme Error XSF Structures (>{extreme_error_percent}% deviation) - Histogram', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)

    # 子图2: 散点图（按 IT 排序，带分隔线）
    extreme_pressures_sorted = [e['total_pressure_gpa'] for e in extreme_errors_sorted_by_it]
    indices = list(range(1, len(extreme_pressures_sorted) + 1))
    colors = ['darkred' if p < extreme_lower or p > extreme_upper else 'orange' for p in extreme_pressures_sorted]

    ax2.axhspan(ax2.get_ylim()[0], extreme_lower, alpha=0.15, color='red', zorder=1)
    ax2.axhspan(extreme_upper, ax2.get_ylim()[1], alpha=0.15, color='red', zorder=1)

    ax2.scatter(indices, extreme_pressures_sorted, c=colors, alpha=0.7, s=30, zorder=3)

    # 绘制 IT 分隔线（蓝色虚线）
    for boundary_idx in it_boundaries:
        ax2.axvline(boundary_idx + 0.5, color='#3498DB', linestyle='--', linewidth=1.5, alpha=0.7, zorder=2)

    # 标注 IT 编号
    if it_boundaries:
        boundaries_with_edges = [0] + it_boundaries + [len(extreme_errors_sorted_by_it)]
        for i in range(len(boundaries_with_edges) - 1):
            start = boundaries_with_edges[i]
            end = boundaries_with_edges[i + 1]
            center = (start + end) / 2 + 1
            it_num = extreme_errors_sorted_by_it[start]['it_number'] if start < len(extreme_errors_sorted_by_it) else -1
            if it_num >= 0:
                ax2.text(center, ax2.get_ylim()[1], f'IT{it_num}',
                        ha='center', va='bottom', fontsize=9, color='#3498DB',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#3498DB', linewidth=1, alpha=0.8))

    ax2.axhline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=2)
    ax2.axhline(extreme_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Upper: {extreme_upper:.2f} GPa', zorder=2)
    ax2.axhline(extreme_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Lower: {extreme_lower:.2f} GPa', zorder=2)

    ax2.set_xlabel('Structure Index', fontsize=12)
    ax2.set_ylabel('Total Pressure (GPa)', fontsize=12)
    ax2.set_title(f'Extreme Error XSF Structures (>{extreme_error_percent}% deviation) - Scatter Plot', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    output_file = 'xsf_extreme_error_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"[INFO] 极大误差结构分布图已保存到: {output_file}")
    print()


def move_extreme_errors_to_trash(extreme_errors_to_move, requested_count, total_count):
    """将极端不合理的 XSF 结构移动到垃圾箱"""
    print("=" * 70)
    print(f"移动极端不合理 XSF 结构到垃圾箱 (前 {requested_count} 个，共 {total_count} 个)")
    print("=" * 70)

    # 创建垃圾箱目录（在当前目录下）
    trash_dir = Path("trashbin")
    trash_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] 垃圾箱目录: {trash_dir.resolve()}")
    print()

    # 显示将要移动的文件
    print(f"[INFO] 将要移动以下 {len(extreme_errors_to_move)} 个结构（按偏差从大到小）:")
    for i, error in enumerate(extreme_errors_to_move, 1):
        xsf_path = error['xsf_path']
        pressure = error['total_pressure_gpa']
        deviation = error['deviation_gpa']
        it_num = error['it_number']
        print(f"  {i:3d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa  |  IT{it_num}")
    print()

    moved_count = 0
    not_found_count = 0

    for error in tqdm(extreme_errors_to_move, desc="移动文件", unit="文件"):
        xsf_path = Path(error['xsf_path'])

        if xsf_path.exists():
            # 构建目标路径，保留原始文件名
            target_path = trash_dir / xsf_path.name

            # 如果目标文件已存在，添加时间戳
            if target_path.exists():
                import time
                timestamp = int(time.time())
                target_path = trash_dir / f"{xsf_path.stem}_{timestamp}{xsf_path.suffix}"

            try:
                shutil.move(str(xsf_path), str(target_path))
                moved_count += 1
            except Exception as e:
                tqdm.write(f"  ✗ 移动失败: {xsf_path}, 错误: {e}")
        else:
            tqdm.write(f"  ✗ 文件不存在: {xsf_path}")
            not_found_count += 1

    print()
    print(f"[INFO] 成功移动: {moved_count} 个文件")
    if not_found_count > 0:
        print(f"[WARNING] 未找到: {not_found_count} 个文件")
    print("=" * 70)
    print()


def read_extreme_error_structures_from_file(file_path='extreme_error_xsf_structures.txt'):
    """从已有的文件中读取极端不合理结构信息"""
    if not Path(file_path).exists():
        print(f"[ERROR] 文件不存在: {file_path}")
        return None

    extreme_errors = []

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                # 解析格式：1. IT9/IT9_SCF_structure.xsf  |  123.45 GPa  |  偏差: +23.45 GPa  |  IT9  [文件不存在]
                match = re.match(r'\s*\d+\.\s+(.+?)\s+\|\s+([-\d.]+)\s+GPa\s+\|\s+偏差:\s+([-+\d.]+)\s+GPa\s+\|\s+IT(\d+)(?:\s+\[文件不存在\])?', line)
                if match:
                    xsf_path = match.group(1)
                    pressure = float(match.group(2))
                    deviation = abs(float(match.group(3)))
                    it_num = int(match.group(4))

                    file_exists = Path(xsf_path).exists()

                    extreme_errors.append({
                        'xsf_path': xsf_path,
                        'total_pressure_gpa': pressure,
                        'deviation_gpa': deviation,
                        'it_number': it_num,
                        'file_exists': file_exists
                    })

        print(f"[INFO] 从 {file_path} 读取到 {len(extreme_errors)} 个极端不合理结构")
        return extreme_errors

    except Exception as e:
        print(f"[ERROR] 读取文件 {file_path} 时出错: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description='检查 XSF 文件中的压强分布，识别并处理偏离目标压强的结构',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 基本分析，目标100 GPa，默认合理阈值±20%
  cd XSF
  python3 check_xsf_pressure.py -p 100

  # 自定义合理阈值为±30%
  python3 check_xsf_pressure.py -p 100 -t 30

  # 过滤极端不合理结构（超过±100%）
  python3 check_xsf_pressure.py -p 100 -e 100

  # 分析并移动偏差最大的前5个极端不合理结构
  python3 check_xsf_pressure.py -p 100 -e 100 -m 5

  # 快速删除模式（读取已有文件）
  python3 check_xsf_pressure.py -r extreme_error_xsf_structures.txt -m 3

输出文件:
  - xsf_pressure_distribution.png: 主压强分布图
  - xsf_extreme_error_distribution.png: 极端不合理结构图
  - outlier_xsf_structures.txt: 偏离合理范围的结构列表
  - extreme_error_xsf_structures.txt: 极端不合理结构列表

注意:
  - XSF 文件必须包含 VIRIAL 矩阵和 PRIMVEC（晶胞矢量）
  - VIRIAL 单位为 eV (电子伏特)，通过 VIRIAL = -stress × volume 关系得到
  - 计算压强：stress = -VIRIAL / volume (单位：eV/Å³)
  - 自动转换为 GPa: 1 eV/Å³ = 160.21766208 GPa
  - 对角元平均值代表静水压强
        """
    )

    parser.add_argument('-p', '--pressure', type=float, default=None,
                        help='目标压强 (GPa) [分析模式必需]')
    parser.add_argument('-t', '--threshold', type=float, default=20.0,
                        help='允许的合理阈值 (百分比%%): 默认 20.0 (±20%%)')
    parser.add_argument('-e', '--extreme-error', type=float, default=None,
                        help='极端不合理阈值 (百分比%%): 超过此阈值的结构单独处理')
    parser.add_argument('-m', '--move-xsf', type=int, default=None, metavar='N',
                        help='移动偏差最大的前 N 个极端不合理结构到 trashbin/')
    parser.add_argument('-r', '--read-file', type=str, default=None, metavar='FILE',
                        help='直接读取已有文件（快速删除模式）')

    args = parser.parse_args()

    print()
    print("=" * 70)
    print("XSF 压强检查工具")
    print("=" * 70)
    print()

    # 判断运行模式
    if args.read_file:
        # 快速删除模式
        print("[INFO] 运行模式: 快速删除模式（读取已有文件）")
        print()

        if args.move_xsf is None:
            parser.error("使用 -r/--read-file 时，必须指定 -m/--move-xsf")

        extreme_errors = read_extreme_error_structures_from_file(args.read_file)

        if extreme_errors is None or len(extreme_errors) == 0:
            print("[ERROR] 无法读取到任何极端不合理结构")
            return

        extreme_errors.sort(key=lambda x: x['deviation_gpa'], reverse=True)

        structures_to_move = extreme_errors[:args.move_xsf]
        move_extreme_errors_to_trash(structures_to_move, args.move_xsf, len(extreme_errors))

        print("=" * 70)
        print("操作完成！")
        print("=" * 70)

    else:
        # 分析模式
        print("[INFO] 运行模式: 分析模式（分析 XSF 文件）")
        print()

        if args.pressure is None:
            parser.error("分析模式下必须指定 -p/--pressure 参数")

        if args.move_xsf is not None and args.extreme_error is None:
            parser.error("--move-xsf 必须与 --extreme-error 一起使用")

        analyze_xsf_pressures(
            target_pressure_gpa=args.pressure,
            threshold_percent=args.threshold,
            extreme_error_percent=args.extreme_error,
            move_extreme_count=args.move_xsf
        )

        print("=" * 70)
        print("分析完成！")
        print("=" * 70)


if __name__ == "__main__":
    main()
