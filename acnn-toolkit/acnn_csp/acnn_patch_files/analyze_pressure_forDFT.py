#!/usr/bin/env python3
"""
压强分析脚本：从DFT计算结果中提取压强信息，分析偏离目标压强的结构

使用方法：在 DFT 目录下运行此脚本
    cd DFT
    python3 analyze_pressure.py -p 100
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


def extract_pressure_from_outcar(outcar_path):
    """
    从OUTCAR文件中提取压强信息
    返回: (external_pressure, pullay_stress) in kB，如果未找到则返回 (None, None)
    """
    external_pressure = None
    pullay_stress = None

    try:
        with open(outcar_path, 'r') as f:
            for line in f:
                # 匹配 "external pressure =       59.36 kB  Pullay stress =     1000.00 kB"
                match = re.search(r'external pressure\s*=\s*([-\d.]+)\s*kB\s+Pullay stress\s*=\s*([-\d.]+)\s*kB', line)
                if match:
                    external_pressure = float(match.group(1))
                    pullay_stress = float(match.group(2))
                    # 取最后一次出现的值（优化完成后的最终值）
    except Exception as e:
        print(f"[WARNING] 无法读取 {outcar_path}: {e}")

    return external_pressure, pullay_stress


def find_corresponding_xsf_path(dft_path):
    """
    根据 IT*/SCF/structure_name/OUTCAR 路径，找到对应的 ../XSF/IT*/structure_name.xsf

    例如：
    IT9/SCF/CeMgH-2x2x2-87345/OUTCAR -> ../XSF/IT9/CeMgH-2x2x2-87345.xsf
    """
    # 解析路径：IT*/SCF/structure_name/OUTCAR
    parts = Path(dft_path).parts

    # 第一个部分应该是 IT*
    if len(parts) < 4 or not parts[0].startswith('IT'):
        return None

    it_name = parts[0]  # IT9
    structure_name = parts[2]  # CeMgH-2x2x2-87345

    # 构建 XSF 路径
    xsf_path = Path("..") / "XSF" / it_name / f"{structure_name}.xsf"

    return xsf_path


def analyze_pressures(target_pressure_gpa, threshold_percent=20.0, extreme_error_percent=None, move_extreme_count=None):
    """
    分析所有OUTCAR文件中的压强信息

    参数:
        target_pressure_gpa: 目标压强 (GPa)
        threshold_percent: 允许的合理阈值 (百分比)
        extreme_error_percent: 极端不合理的最小阈值 (百分比)，超过此阈值的结构将被单独处理
        move_extreme_count: 移动极端不合理结构的前 n 个到垃圾箱（None 表示不移动）
    """
    # 检查当前目录
    cwd = Path.cwd()
    if cwd.name != "DFT":
        print(f"[ERROR] 当前目录不是 DFT 目录")
        print(f"[ERROR] 当前目录: {cwd}")
        print(f"[ERROR] 请先切换到 DFT 目录: cd DFT")
        return

    # 计算合理阈值范围（基于百分比）
    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    print(f"[INFO] 当前工作目录: {cwd}")
    print(f"[INFO] 目标压强: {target_pressure_gpa} GPa")
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

    # 查找所有 OUTCAR 文件（相对于 DFT 目录）
    outcar_pattern = "IT*/SCF/*/OUTCAR"
    outcar_files = glob.glob(outcar_pattern, recursive=False)

    if not outcar_files:
        print(f"[ERROR] 未找到任何匹配的 OUTCAR 文件: {outcar_pattern}")
        print(f"[ERROR] 请确认：")
        print(f"[ERROR]   1. 当前在 DFT 目录下")
        print(f"[ERROR]   2. 存在 IT*/SCF/*/OUTCAR 文件")
        return

    print(f"[INFO] 找到 {len(outcar_files)} 个 OUTCAR 文件")
    print()

    # 存储所有压强数据
    pressure_data = []  # 正常数据（在极大误差范围内，用于主图）
    outliers = []  # 偏离目标的结构（超过threshold但在extreme范围内）
    extreme_errors = []  # 极大误差的结构（超过extreme范围）

    # 提取压强信息
    print("[INFO] 正在提取 OUTCAR 文件中的压强信息...")
    for outcar_path in tqdm(outcar_files, desc="处理 OUTCAR", unit="文件"):
        ext_press, pullay = extract_pressure_from_outcar(outcar_path)

        if ext_press is not None and pullay is not None:
            total_pressure_kb = ext_press + pullay
            total_pressure_gpa = total_pressure_kb / 10.0

            data_entry = {
                'outcar_path': outcar_path,
                'external_pressure': ext_press,
                'pullay_stress': pullay,
                'total_pressure_kb': total_pressure_kb,
                'total_pressure_gpa': total_pressure_gpa
            }

            # 检查是否为极大误差
            is_extreme = False
            if extreme_error_percent is not None:
                extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
                extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)
                if total_pressure_gpa < extreme_lower or total_pressure_gpa > extreme_upper:
                    is_extreme = True
                    xsf_path = find_corresponding_xsf_path(outcar_path)
                    extreme_errors.append({
                        'outcar_path': outcar_path,
                        'xsf_path': xsf_path,
                        'total_pressure_gpa': total_pressure_gpa,
                        'deviation_gpa': abs(total_pressure_gpa - target_pressure_gpa)
                    })

            # 如果不是极大误差，加入主数据
            if not is_extreme:
                pressure_data.append(data_entry)

                # 检查是否偏离合理范围
                if total_pressure_gpa < threshold_lower or total_pressure_gpa > threshold_upper:
                    xsf_path = find_corresponding_xsf_path(outcar_path)
                    outliers.append({
                        'outcar_path': outcar_path,
                        'xsf_path': xsf_path,
                        'total_pressure_gpa': total_pressure_gpa,
                        'deviation_gpa': abs(total_pressure_gpa - target_pressure_gpa)
                    })

    if not pressure_data and not extreme_errors:
        print("[ERROR] 未能从任何 OUTCAR 文件中提取到压强信息")
        return

    print(f"[INFO] 成功提取 {len(pressure_data) + len(extreme_errors)} 个结构的压强信息")
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
        print("=" * 70)
        print(f"偏离合理范围（±{threshold_percent}%）的结构列表")
        print("=" * 70)

        # 按偏差从大到小排序
        outliers.sort(key=lambda x: x['deviation_gpa'], reverse=True)

        # 输出到文件
        with open('outlier_structures.txt', 'w') as f:
            f.write(f"# 偏离合理范围的结构 (目标压强 {target_pressure_gpa} GPa, 合理阈值 ±{threshold_percent}%)\n")
            f.write(f"# 共 {len(outliers)} 个结构\n")
            f.write(f"# 合理范围: {threshold_lower:.2f} ~ {threshold_upper:.2f} GPa\n")
            f.write("#\n")
            f.write("# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa)\n")
            f.write("#" + "=" * 68 + "\n\n")

            for i, outlier in enumerate(outliers, 1):
                xsf_path = outlier['xsf_path']
                pressure = outlier['total_pressure_gpa']
                deviation = outlier['deviation_gpa']

                line = f"{i:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa\n"
                print(f"  {line.strip()}")
                f.write(line)

        print()
        print(f"[INFO] 偏离结构列表已保存到: outlier_structures.txt")
        print()

    # 处理极大误差的结构
    if extreme_errors:
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

        # 输出到文件
        with open('extreme_error_structures.txt', 'w') as f:
            f.write(f"# 极大误差结构列表 (超过目标压强 {target_pressure_gpa} GPa 的 ±{extreme_error_percent}%)\n")
            f.write(f"# 共 {len(extreme_errors)} 个结构\n")
            f.write(f"# 极大误差范围: < {target_pressure_gpa * (1 - extreme_error_percent / 100.0):.2f} GPa 或 > {target_pressure_gpa * (1 + extreme_error_percent / 100.0):.2f} GPa\n")
            f.write("#\n")
            f.write("# 格式: XSF路径 | 实际压强(GPa) | 偏差(GPa)\n")
            f.write("#" + "=" * 68 + "\n\n")

            for i, error in enumerate(extreme_errors, 1):
                xsf_path = error['xsf_path']
                pressure = error['total_pressure_gpa']
                deviation = error['deviation_gpa']

                line = f"{i:4d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa\n"
                if i <= 10:  # 只打印前10个到终端
                    print(f"  {line.strip()}")
                f.write(line)

            if len(extreme_errors) > 10:
                print(f"  ... (还有 {len(extreme_errors) - 10} 个，详见文件)")

        print()
        print(f"[INFO] 极大误差结构列表已保存到: extreme_error_structures.txt")
        print()

        # 如果指定了移动参数，移动前 n 个极端不合理的结构
        if move_extreme_count is not None and move_extreme_count > 0:
            structures_to_move = extreme_errors[:move_extreme_count]
            move_extreme_errors_to_trash(structures_to_move, move_extreme_count, len(extreme_errors))

        # 绘制极大误差结构的分布图
        plot_extreme_error_distribution(extreme_errors, target_pressure_gpa, extreme_error_percent)


def plot_pressure_distribution(pressure_data, target_pressure_gpa, threshold_percent):
    """
    绘制压强分布图

    参数:
        pressure_data: 压强数据列表
        target_pressure_gpa: 目标压强 (GPa)
        threshold_percent: 合理阈值 (百分比)
    """
    total_pressures = [d['total_pressure_gpa'] for d in pressure_data]

    # 计算合理范围
    threshold_lower = target_pressure_gpa * (1 - threshold_percent / 100.0)
    threshold_upper = target_pressure_gpa * (1 + threshold_percent / 100.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # 子图1: 直方图
    n, bins, patches = ax1.hist(total_pressures, bins=50, alpha=0.7, color='steelblue', edgecolor='black')

    # 填充合理范围（先画，作为背景）
    ax1.axvspan(threshold_lower, threshold_upper,
                alpha=0.15, color='green', label=f'Reasonable Range (±{threshold_percent}%)', zorder=1)

    # 标记目标压强（红色粗线）
    ax1.axvline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=3)

    # 标记阈值范围上限（橙色粗虚线）
    ax1.axvline(threshold_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Upper Limit: {threshold_upper:.2f} GPa', zorder=2)

    # 标记阈值范围下限（橙色粗虚线）
    ax1.axvline(threshold_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Lower Limit: {threshold_lower:.2f} GPa', zorder=2)

    # 在图上添加文本标注
    y_max = n.max() * 0.95
    ax1.text(target_pressure_gpa, y_max, f'{target_pressure_gpa} GPa',
             ha='center', va='top', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))
    ax1.text(threshold_lower, y_max * 0.75, f'{threshold_lower:.2f} GPa',
             ha='center', va='top', fontsize=9, color='#E67E22',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E67E22', linewidth=1.5))
    ax1.text(threshold_upper, y_max * 0.75, f'{threshold_upper:.2f} GPa',
             ha='center', va='top', fontsize=9, color='#E67E22',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E67E22', linewidth=1.5))

    ax1.set_xlabel('Total Pressure (GPa)', fontsize=12)
    ax1.set_ylabel('Number of Structures', fontsize=12)
    ax1.set_title('Pressure Distribution of DFT Structures (Histogram)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)

    # 子图2: 散点图（按序号排列）
    indices = list(range(1, len(total_pressures) + 1))
    colors = ['green' if threshold_lower <= p <= threshold_upper else 'red' for p in total_pressures]

    # 填充合理范围（先画，作为背景）
    ax2.axhspan(threshold_lower, threshold_upper,
                alpha=0.15, color='green', label=f'Reasonable Range (±{threshold_percent}%)', zorder=1)

    # 绘制散点
    ax2.scatter(indices, total_pressures, c=colors, alpha=0.6, s=20, zorder=3)

    # 标记目标压强（红色粗线）
    ax2.axhline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=2)

    # 标记阈值范围上限（橙色粗虚线）
    ax2.axhline(threshold_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Upper Limit: {threshold_upper:.2f} GPa', zorder=2)

    # 标记阈值范围下限（橙色粗虚线）
    ax2.axhline(threshold_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Lower Limit: {threshold_lower:.2f} GPa', zorder=2)

    # 在右侧添加文本标注
    x_max = len(indices) * 0.98
    ax2.text(x_max, target_pressure_gpa, f'{target_pressure_gpa} GPa',
             ha='right', va='center', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))
    ax2.text(x_max, threshold_lower, f'{threshold_lower:.2f} GPa',
             ha='right', va='center', fontsize=9, color='#E67E22',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E67E22', linewidth=1.5))
    ax2.text(x_max, threshold_upper, f'{threshold_upper:.2f} GPa',
             ha='right', va='center', fontsize=9, color='#E67E22',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E67E22', linewidth=1.5))

    ax2.set_xlabel('Structure Index', fontsize=12)
    ax2.set_ylabel('Total Pressure (GPa)', fontsize=12)
    ax2.set_title('Pressure Distribution of DFT Structures (Scatter Plot)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # 保存图片
    output_file = 'pressure_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"[INFO] 压强分布图已保存到: {output_file}")
    print()

    # 如果在非交互环境，关闭图形
    # plt.close()


def plot_extreme_error_distribution(extreme_errors, target_pressure_gpa, extreme_error_percent):
    """
    绘制极大误差结构的压强分布图
    """
    extreme_pressures = [e['total_pressure_gpa'] for e in extreme_errors]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # 计算极大误差范围
    extreme_lower = target_pressure_gpa * (1 - extreme_error_percent / 100.0)
    extreme_upper = target_pressure_gpa * (1 + extreme_error_percent / 100.0)

    # 子图1: 直方图
    n, bins, patches = ax1.hist(extreme_pressures, bins=30, alpha=0.7, color='crimson', edgecolor='black')

    # 标记目标压强（红色粗线）
    ax1.axvline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=3)

    # 标记极大误差范围边界（橙色粗虚线）
    ax1.axvline(extreme_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Upper: {extreme_upper:.2f} GPa', zorder=2)
    ax1.axvline(extreme_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Lower: {extreme_lower:.2f} GPa', zorder=2)

    # 填充极大误差区域（浅红色）
    ax1.axvspan(ax1.get_xlim()[0], extreme_lower, alpha=0.15, color='red', label='Extreme Error Zone', zorder=1)
    ax1.axvspan(extreme_upper, ax1.get_xlim()[1], alpha=0.15, color='red', zorder=1)

    # 在图上添加文本标注
    y_max = n.max() * 0.95
    ax1.text(target_pressure_gpa, y_max, f'{target_pressure_gpa} GPa',
             ha='center', va='top', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))

    ax1.set_xlabel('Total Pressure (GPa)', fontsize=12)
    ax1.set_ylabel('Number of Structures', fontsize=12)
    ax1.set_title(f'Extreme Error Structures (>{extreme_error_percent}% deviation) - Histogram', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='best')
    ax1.grid(True, alpha=0.3)

    # 子图2: 散点图（按序号排列）
    indices = list(range(1, len(extreme_pressures) + 1))
    colors = ['darkred' if p < extreme_lower or p > extreme_upper else 'orange' for p in extreme_pressures]

    # 填充极大误差区域（浅红色）
    ax2.axhspan(ax2.get_ylim()[0], extreme_lower, alpha=0.15, color='red', zorder=1)
    ax2.axhspan(extreme_upper, ax2.get_ylim()[1], alpha=0.15, color='red', zorder=1)

    # 绘制散点
    ax2.scatter(indices, extreme_pressures, c=colors, alpha=0.7, s=30, zorder=3)

    # 标记目标压强（红色粗线）
    ax2.axhline(target_pressure_gpa, color='#E74C3C', linestyle='-', linewidth=3,
                label=f'Target: {target_pressure_gpa} GPa', zorder=2)

    # 标记极大误差范围边界（橙色粗虚线）
    ax2.axhline(extreme_upper, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Upper: {extreme_upper:.2f} GPa', zorder=2)
    ax2.axhline(extreme_lower, color='#E67E22', linestyle='--', linewidth=2.5,
                label=f'Extreme Lower: {extreme_lower:.2f} GPa', zorder=2)

    # 在右侧添加文本标注
    x_max = len(indices) * 0.98
    ax2.text(x_max, target_pressure_gpa, f'{target_pressure_gpa} GPa',
             ha='right', va='center', fontsize=11, fontweight='bold',
             color='#E74C3C', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C', linewidth=2))

    ax2.set_xlabel('Structure Index', fontsize=12)
    ax2.set_ylabel('Total Pressure (GPa)', fontsize=12)
    ax2.set_title(f'Extreme Error Structures (>{extreme_error_percent}% deviation) - Scatter Plot', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # 保存图片
    output_file = 'extreme_error_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"[INFO] 极大误差结构分布图已保存到: {output_file}")
    print()

    # 如果在非交互环境，关闭图形
    # plt.close()


def move_extreme_errors_to_trash(extreme_errors_to_move, requested_count, total_count):
    """
    将极端不合理的结构移动到垃圾箱

    参数:
        extreme_errors_to_move: 要移动的极端不合理结构列表（已排序的前 n 个）
        requested_count: 用户请求移动的数量
        total_count: 极端不合理结构总数
    """
    print("=" * 70)
    print(f"移动极端不合理结构到垃圾箱 (前 {requested_count} 个，共 {total_count} 个)")
    print("=" * 70)

    # 创建垃圾箱目录（在 DFT 的同级目录 XSF 下）
    trash_dir = Path("../XSF/trashbin")
    trash_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] 垃圾箱目录: {trash_dir.resolve()}")
    print()

    # 显示将要移动的文件
    print(f"[INFO] 将要移动以下 {len(extreme_errors_to_move)} 个结构（按偏差从大到小）:")
    for i, error in enumerate(extreme_errors_to_move, 1):
        xsf_path = error['xsf_path']
        pressure = error['total_pressure_gpa']
        deviation = error['deviation_gpa']
        print(f"  {i:3d}. {xsf_path}  |  {pressure:7.2f} GPa  |  偏差: {deviation:+7.2f} GPa")
    print()

    moved_count = 0
    not_found_count = 0

    for error in tqdm(extreme_errors_to_move, desc="移动文件", unit="文件"):
        xsf_path = error['xsf_path']

        if xsf_path and Path(xsf_path).exists():
            # 构建目标路径
            # 保留 IT* 信息，例如: ../XSF/trashbin/IT9_structure.xsf
            xsf_file = Path(xsf_path)
            it_name = xsf_file.parent.name  # IT9
            new_name = f"{it_name}_{xsf_file.name}"
            target_path = trash_dir / new_name

            # 移动文件
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


def main():
    parser = argparse.ArgumentParser(
        description='分析DFT计算结构的压强分布，识别并处理偏离目标压强的结构',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 基本分析，目标100 GPa，默认合理阈值±20%
  cd DFT
  python3 analyze_pressure.py -p 100

  # 自定义合理阈值为±30%（70~130 GPa为合理范围）
  python3 analyze_pressure.py -p 100 -t 30

  # 过滤极端不合理结构（超过±100%的结构单独处理）
  python3 analyze_pressure.py -p 100 --extreme-error 100
  # <0 或 >200 GPa 的结构不在主图中显示

  # 组合使用：合理阈值±30%，极端阈值±100%
  python3 analyze_pressure.py -p 100 -t 30 --extreme-error 100
  # 70~130 GPa: 合理范围（绿色）
  # 0~70 和 130~200 GPa: 偏离但在主图中（红色）
  # <0 或 >200 GPa: 极端不合理，单独图表

  # 分析并将偏差最大的前5个极端不合理结构移动到垃圾箱
  python3 analyze_pressure.py -p 100 --extreme-error 100 --move-xsf 5

输出文件:
  - pressure_distribution.png: 主压强分布图
  - extreme_error_distribution.png: 极端不合理结构图（如果有）
  - outlier_structures.txt: 主图中偏离合理范围的结构列表
  - extreme_error_structures.txt: 极端不合理结构列表（如果有）

注意:
  - --move-xsf 必须与 --extreme-error 一起使用
  - 只会移动极端不合理的结构（按偏差从大到小排序后的前 n 个）
  - 移动前会显示文件列表供确认
        """
    )

    parser.add_argument('-p', '--pressure', type=float, required=True,
                        help='目标压强 (GPa)')
    parser.add_argument('-t', '--threshold', type=float, default=20.0,
                        help='允许的合理阈值 (百分比%%): 在此范围内的结构视为合理。默认: 20.0 (表示±20%%)')
    parser.add_argument('--extreme-error', type=float, default=None,
                        help='极端不合理的最小阈值 (百分比%%): 超过此阈值的结构不在主图中显示，单独处理。例如: 100 表示超过±100%%的结构是极端不合理的')
    parser.add_argument('--move-xsf', type=int, default=None, metavar='N',
                        help='移动偏差最大的前 N 个极端不合理结构从 ../XSF/IT* 到 ../XSF/trashbin (必须与 --extreme-error 一起使用)')

    args = parser.parse_args()

    # 验证参数
    if args.move_xsf is not None and args.extreme_error is None:
        parser.error("--move-xsf 必须与 --extreme-error 一起使用")

    print()
    print("=" * 70)
    print("DFT 压强分析工具")
    print("=" * 70)
    print()

    analyze_pressures(
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
