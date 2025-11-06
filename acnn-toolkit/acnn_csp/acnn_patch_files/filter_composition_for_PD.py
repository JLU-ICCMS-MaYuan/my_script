#!/usr/bin/env python3
"""
配比筛选工具：统计、筛选、移动不满足配比约束条件的结构

功能：
1. 统计配比分布（保留原有功能）
2. 根据元素比例约束筛选结构
3. 批量移动不满足条件的结构到垃圾箱

使用示例：
    # 统计配比分布
    python filter_composition.py IT0 --stats

    # 筛选满足 Ce:Mg <= 5:1 的结构
    python filter_composition.py IT0 --ratio-limit Ce:Mg:5 --filter

    # 移动所有IT*中不满足条件的结构到垃圾箱
    python filter_composition.py --ratio-limit Ce:Mg:5 --move-invalid

    # 模拟运行（不实际移动）
    python filter_composition.py --ratio-limit Ce:Mg:5 --move-invalid --dry-run
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from collections import Counter, OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def parse_iteration_label(label: str) -> str:
    """解析迭代编号，统一转换为 IT* 格式"""
    value = label.strip()
    if not value:
        raise ValueError("迭代编号不能为空")
    upper = value.upper()
    if upper.startswith("IT"):
        return upper
    return f"IT{upper}"


def gcd_many(values: Iterable[int]) -> int:
    """计算多个整数的最大公约数"""
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration:
        return 1
    result = abs(result)
    for value in iterator:
        result = math.gcd(result, abs(value))
        if result == 1:
            return 1
    return max(result, 1)


def parse_composition(path: Path) -> OrderedDict[str, int]:
    """
    从 .res 文件解析元素配比

    返回：OrderedDict，按 SFAC 顺序排列，例如 {"Ce": 1, "Mg": 10, "H": 22}
    """
    counts: Counter[str] = Counter()
    sfac_order: List[str] = []
    in_atoms = False

    try:
        with path.open("r", encoding="ascii", errors="ignore") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                upper = line.upper()

                # 跳过标题行
                if upper.startswith("TITL") or upper.startswith("CELL") or upper.startswith("LATT"):
                    continue

                # 读取 SFAC（元素顺序）
                if upper.startswith("SFAC"):
                    parts = raw_line.split()
                    for token in parts[1:]:
                        token_clean = token.strip()
                        if not token_clean:
                            continue
                        sfac_order.append(token_clean)
                    in_atoms = True
                    continue

                # 跳过其他元数据
                if upper.startswith("UNIT") or upper.startswith("REM"):
                    continue

                # 结束标记
                if upper == "END":
                    break

                # 解析原子数据
                if not in_atoms:
                    continue

                parts = raw_line.split()
                if not parts:
                    continue

                element = parts[0].strip()
                if not element or element.upper() == "END":
                    continue

                counts[element] += 1
    except Exception as e:
        print(f"警告：解析文件 {path} 时出错: {e}", file=sys.stderr)
        return OrderedDict()

    # 按 SFAC 顺序构建 OrderedDict
    ordered_counts: OrderedDict[str, int] = OrderedDict()
    if counts:
        if not sfac_order:
            # 若未能读取到 SFAC 顺序，则按元素名称排序
            for elem in sorted(counts.keys()):
                ordered_counts[elem] = counts[elem]
        else:
            added = set()
            for elem in sfac_order:
                if elem in counts and elem not in added:
                    ordered_counts[elem] = counts[elem]
                    added.add(elem)
            # 补充在 SFAC 中未出现但统计到的元素
            for elem in counts:
                if elem not in ordered_counts:
                    ordered_counts[elem] = counts[elem]

    return ordered_counts


def composition_signature(ordered_counts: OrderedDict[str, int]) -> str:
    """
    生成配比签名（化约到最简整数比）

    例如：{"Ce": 2, "Mg": 10, "H": 22} -> "CeMg5H11"
    """
    if not ordered_counts:
        return "未知配比"

    gcd_value = gcd_many(ordered_counts.values())
    parts: List[str] = []
    for elem, count in ordered_counts.items():
        base = count // gcd_value
        if base == 1:
            parts.append(elem)
        else:
            parts.append(f"{elem}{base}")

    return "".join(parts)


def parse_ratio_limit(ratio_str: str) -> Tuple[str, str, float]:
    """
    解析比例约束参数

    格式：ELEM1:ELEM2:MAX_RATIO
    例如："Ce:Mg:5" -> ("Ce", "Mg", 5.0)
    """
    parts = ratio_str.split(":")
    if len(parts) != 3:
        raise ValueError(f"比例约束格式错误：'{ratio_str}'，应为 ELEM1:ELEM2:MAX_RATIO")

    elem1, elem2, ratio_str = parts
    elem1 = elem1.strip()
    elem2 = elem2.strip()

    if not elem1 or not elem2:
        raise ValueError(f"元素名称不能为空：'{ratio_str}'")

    try:
        max_ratio = float(ratio_str.strip())
    except ValueError:
        raise ValueError(f"比例值必须是数字：'{ratio_str}'")

    if max_ratio <= 0:
        raise ValueError(f"比例值必须大于0：'{ratio_str}'")

    return elem1, elem2, max_ratio


def check_ratio_constraint(
    composition: OrderedDict[str, int],
    elem1: str,
    elem2: str,
    max_ratio: float
) -> Tuple[bool, Optional[float]]:
    """
    检查两个元素的比例是否满足约束

    约束规则：max(count1, count2) / min(count1, count2) <= max_ratio

    返回：(是否满足, 实际比例)
    """
    count1 = composition.get(elem1, 0)
    count2 = composition.get(elem2, 0)

    # 如果某个元素不存在，视为不满足条件
    if count1 == 0 or count2 == 0:
        return False, None

    # 计算比例：大的除以小的
    actual_ratio = max(count1, count2) / min(count1, count2)

    # 判断是否满足约束
    is_valid = actual_ratio <= max_ratio

    return is_valid, actual_ratio


def format_ratio(composition: OrderedDict[str, int], elem1: str, elem2: str) -> str:
    """格式化显示比例，例如 "Ce:Mg = 1:10" """
    count1 = composition.get(elem1, 0)
    count2 = composition.get(elem2, 0)

    if count1 == 0 or count2 == 0:
        return f"{elem1}:{elem2} = N/A"

    # 化简比例
    gcd_val = math.gcd(count1, count2)
    simplified1 = count1 // gcd_val
    simplified2 = count2 // gcd_val

    return f"{elem1}:{elem2} = {simplified1}:{simplified2}"


def get_element_count(composition: OrderedDict[str, int]) -> int:
    """获取配比中的元素数量"""
    return len([elem for elem, count in composition.items() if count > 0])


def should_exclude_by_element_count(
    composition: OrderedDict[str, int],
    exclude_counts: List[int]
) -> bool:
    """
    根据元素数量判断是否应该排除

    例如：exclude_counts = [1, 2, 3] 表示单质、二元、三元化合物不移动
    """
    if not exclude_counts:
        return False

    elem_count = get_element_count(composition)
    return elem_count in exclude_counts


def check_all_constraints(
    composition: OrderedDict[str, int],
    constraints: List[Tuple[str, str, float]],
    exclude_element_counts: Optional[List[int]] = None
) -> Tuple[bool, List[str]]:
    """
    检查所有比例约束

    返回：(是否全部满足, 不满足的约束描述列表)
    """
    # 先检查是否因元素数量而排除
    if exclude_element_counts and should_exclude_by_element_count(composition, exclude_element_counts):
        # 被排除的配比视为"满足条件"（不移动）
        return True, []

    all_valid = True
    violations: List[str] = []

    for elem1, elem2, max_ratio in constraints:
        is_valid, actual_ratio = check_ratio_constraint(composition, elem1, elem2, max_ratio)

        if not is_valid:
            all_valid = False
            ratio_str = format_ratio(composition, elem1, elem2)
            if actual_ratio is not None:
                violations.append(f"{ratio_str} > {max_ratio}:1")
            else:
                violations.append(f"{elem1}或{elem2}不存在")

    return all_valid, violations


def find_it_directories(pd_dir: Path, iteration: Optional[str] = None) -> List[Path]:
    """
    查找 IT* 目录

    如果指定 iteration，只返回该目录；否则返回所有 IT* 目录
    """
    if iteration:
        it_dir = pd_dir / iteration
        if not it_dir.is_dir():
            return []
        return [it_dir]
    else:
        # 查找所有 IT* 目录
        it_dirs = sorted([d for d in pd_dir.glob("IT*") if d.is_dir()])
        return it_dirs


def show_progress_bar(current: int, total: int, bar_width: int = 40):
    """显示进度条"""
    if not sys.stdout.isatty():
        return

    filled = int(bar_width * current / total) if total > 0 else 0
    bar = "#" * filled + "." * (bar_width - filled)
    sys.stdout.write(f"\r解析进度：[{bar}] {current}/{total}")
    sys.stdout.flush()


def move_single_file(args_tuple) -> Tuple[bool, Path, Path, str]:
    """
    移动单个文件（用于并行处理）

    参数：(源路径, 垃圾箱目录, IT名称)
    返回：(是否成功, 源路径, 目标路径或错误信息, 原因)
    """
    src_path, trash_dir, it_name, reason = args_tuple

    # 构建目标路径，保留 IT 目录结构
    target_dir = trash_dir / it_name
    target_path = target_dir / src_path.name

    try:
        # 创建目标目录
        target_dir.mkdir(parents=True, exist_ok=True)

        # 移动文件（直接覆盖如果存在）
        shutil.move(str(src_path), str(target_path))

        return True, src_path, target_path, reason
    except Exception as e:
        return False, src_path, str(e), reason


def restore_single_file(args_tuple) -> Tuple[bool, Path, Path, str]:
    """
    还原单个文件（用于并行处理）

    参数：(垃圾箱中的路径, 原始路径)
    返回：(是否成功, 垃圾箱路径, 原始路径或错误信息, "")
    """
    trash_path, original_path = args_tuple

    try:
        # 确保原始目录存在
        original_path.parent.mkdir(parents=True, exist_ok=True)

        # 移动文件回原位置（直接覆盖如果存在）
        shutil.move(str(trash_path), str(original_path))

        return True, trash_path, original_path, ""
    except Exception as e:
        return False, trash_path, str(e), ""


def write_move_log(
    log_path: Path,
    moved_files: List[Tuple[Path, Path, str]],
    constraints: List[Tuple[str, str, float]]
):
    """
    写入移动操作日志

    moved_files: [(源路径, 目标路径, 违反原因), ...]
    """
    with log_path.open("w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("配比筛选 - 文件移动日志\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"移动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("筛选条件：\n")
        for elem1, elem2, max_ratio in constraints:
            f.write(f"  - {elem1}:{elem2} 比例 <= {max_ratio}:1\n")
        f.write("\n")

        f.write(f"移动文件数量：{len(moved_files)}\n\n")

        f.write("移动详情：\n")
        f.write("-" * 70 + "\n")
        for src, dst, reason in moved_files:
            f.write(f"源文件：{src}\n")
            f.write(f"目标：  {dst}\n")
            f.write(f"原因：  {reason}\n")
            f.write("-" * 70 + "\n")


def analyze_element_ranges(
    all_compositions: Dict[Path, OrderedDict[str, int]]
) -> Dict[str, Dict]:
    """
    分析每个元素的数量范围和出现的所有可能值

    返回：{
        'Ce': {
            'min': 1, 'max': 10, 'avg': 2.3,
            'values': {1, 2, 3, 5, 10},
            'count': 150
        },
        ...
    }
    """
    element_data: Dict[str, List[int]] = {}

    # 收集所有元素的数量
    for composition in all_compositions.values():
        for elem, count in composition.items():
            if elem not in element_data:
                element_data[elem] = []
            element_data[elem].append(count)

    # 计算统计信息
    element_stats: Dict[str, Dict] = {}
    for elem, counts in element_data.items():
        element_stats[elem] = {
            'min': min(counts),
            'max': max(counts),
            'avg': sum(counts) / len(counts),
            'values': sorted(set(counts)),
            'count': len(counts)
        }

    return element_stats


def analyze_element_combinations(
    all_compositions: Dict[Path, OrderedDict[str, int]]
) -> Dict[str, Dict[int, List[int]]]:
    """
    分析元素组合的覆盖情况

    返回：{
        'Ce-Mg': {
            1: [5, 8, 10],      # Ce=1 时，Mg 的所有值
            2: [8, 10, 15],     # Ce=2 时，Mg 的所有值
            ...
        },
        ...
    }
    """
    # 获取所有出现的元素
    all_elements = set()
    for composition in all_compositions.values():
        all_elements.update(composition.keys())

    all_elements = sorted(all_elements)

    # 分析两两元素的组合
    combinations: Dict[str, Dict[int, set]] = {}

    if len(all_elements) >= 2:
        # 只分析前两个主要元素（通常是金属元素）
        for i in range(len(all_elements)):
            for j in range(i + 1, len(all_elements)):
                elem1, elem2 = all_elements[i], all_elements[j]
                key = f"{elem1}-{elem2}"
                combinations[key] = {}

                for composition in all_compositions.values():
                    count1 = composition.get(elem1, 0)
                    count2 = composition.get(elem2, 0)

                    if count1 > 0 and count2 > 0:
                        if count1 not in combinations[key]:
                            combinations[key][count1] = set()
                        combinations[key][count1].add(count2)

    # 转换 set 为 sorted list
    result: Dict[str, Dict[int, List[int]]] = {}
    for key, data in combinations.items():
        result[key] = {k: sorted(v) for k, v in data.items()}

    return result


def generate_detailed_report(
    it_dirs: List[Path],
    total_files: int,
    all_compositions: Dict[Path, OrderedDict[str, int]],
    composition_counter: Counter[str],
    constraints: List[Tuple[str, str, float]],
    valid_counter: Optional[Counter[str]] = None,
    invalid_counter: Optional[Counter[str]] = None
) -> str:
    """生成详细的配比分析报告"""
    lines: List[str] = []

    # 标题
    lines.append("=" * 70)
    lines.append("配比统计详细报告")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"扫描目录：{', '.join(d.name for d in it_dirs)}")
    lines.append(f"文件总数：{total_files}")
    lines.append("")

    # 筛选条件
    if constraints:
        lines.append("筛选条件：")
        for elem1, elem2, max_ratio in constraints:
            lines.append(f"  - {elem1}:{elem2} 比例 <= {max_ratio}:1")
        lines.append("")
        valid_count = sum(valid_counter.values()) if valid_counter else 0
        invalid_count = sum(invalid_counter.values()) if invalid_counter else 0
        lines.append(f"满足条件：{valid_count} ({valid_count*100/total_files:.1f}%)")
        lines.append(f"不满足条件：{invalid_count} ({invalid_count*100/total_files:.1f}%)")
        lines.append("")

    # 配比分布
    lines.append("=" * 70)
    lines.append("所有配比分布（按数量排序）")
    lines.append("=" * 70)
    lines.append(f"{'配比':<25} {'数量':>10}  {'百分比':>8}")
    lines.append("-" * 70)

    for signature, count in composition_counter.most_common():
        percentage = count * 100 / total_files
        lines.append(f"{signature:<25} {count:>10}  {percentage:>7.2f}%")

    lines.append("")

    # 元素数量范围统计
    lines.append("=" * 70)
    lines.append("元素数量范围统计")
    lines.append("=" * 70)

    element_stats = analyze_element_ranges(all_compositions)

    lines.append(f"{'元素':<8} {'最小值':>8} {'最大值':>8} {'平均值':>10} {'出现次数':>10}")
    lines.append("-" * 70)

    for elem in sorted(element_stats.keys()):
        stats = element_stats[elem]
        lines.append(
            f"{elem:<8} {stats['min']:>8} {stats['max']:>8} "
            f"{stats['avg']:>10.2f} {stats['count']:>10}"
        )

    lines.append("")

    # 每个元素的所有出现值
    lines.append("=" * 70)
    lines.append("每个元素的所有可能值")
    lines.append("=" * 70)

    for elem in sorted(element_stats.keys()):
        values = element_stats[elem]['values']
        lines.append(f"{elem}: {', '.join(map(str, values))}")

    lines.append("")

    # 元素组合分析
    lines.append("=" * 70)
    lines.append("元素组合分析（已探索的配比组合）")
    lines.append("=" * 70)

    combinations = analyze_element_combinations(all_compositions)

    for combo_key in sorted(combinations.keys()):
        combo_data = combinations[combo_key]
        elem1, elem2 = combo_key.split('-')

        lines.append(f"\n{combo_key} 组合:")
        lines.append("-" * 70)

        for val1 in sorted(combo_data.keys()):
            val2_list = combo_data[val1]
            lines.append(f"  {elem1}={val1}: {elem2}={', '.join(map(str, val2_list))}")

    lines.append("")

    # 配比空隙分析（未探索的组合提示）
    lines.append("=" * 70)
    lines.append("配比空隙分析（可能未探索的区域）")
    lines.append("=" * 70)

    for combo_key in sorted(combinations.keys()):
        combo_data = combinations[combo_key]
        elem1, elem2 = combo_key.split('-')

        elem1_values = sorted(combo_data.keys())
        if not elem1_values:
            continue

        lines.append(f"\n{combo_key} 可能的空隙:")
        lines.append("-" * 70)

        # 检查连续性
        missing_combos: List[Tuple[int, int]] = []

        for val1 in elem1_values:
            val2_list = combo_data[val1]
            val2_min, val2_max = min(val2_list), max(val2_list)

            # 检查 val2 范围内的空隙
            for val2 in range(val2_min, val2_max + 1):
                if val2 not in val2_list:
                    missing_combos.append((val1, val2))

        if missing_combos:
            for val1, val2 in missing_combos[:20]:  # 只显示前20个
                lines.append(f"  {elem1}={val1}, {elem2}={val2} (未发现)")
            if len(missing_combos) > 20:
                lines.append(f"  ... (还有 {len(missing_combos) - 20} 个组合未列出)")
        else:
            lines.append(f"  在已有范围内未发现明显空隙")

    lines.append("")

    # 满足/不满足条件的配比分布（如果有筛选条件）
    if constraints and valid_counter and invalid_counter:
        lines.append("=" * 70)
        lines.append("满足条件的配比分布")
        lines.append("=" * 70)
        lines.append(f"{'配比':<25} {'数量':>10}")
        lines.append("-" * 70)

        for signature, count in valid_counter.most_common():
            lines.append(f"{signature:<25} {count:>10}")

        lines.append("")
        lines.append("=" * 70)
        lines.append("不满足条件的配比分布")
        lines.append("=" * 70)
        lines.append(f"{'配比':<25} {'数量':>10}")
        lines.append("-" * 70)

        for signature, count in invalid_counter.most_common():
            lines.append(f"{signature:<25} {count:>10}")

        lines.append("")

    lines.append("=" * 70)
    lines.append("报告结束")
    lines.append("=" * 70)

    return "\n".join(lines)


def stats_mode(
    it_dirs: List[Path],
    constraints: List[Tuple[str, str, float]],
    output_file: Optional[str],
    exclude_element_counts: Optional[List[int]] = None
):
    """统计模式：显示配比分布"""
    print("\n" + "=" * 70)
    print("配比统计模式")
    print("=" * 70 + "\n")

    if constraints:
        print("筛选条件：")
        for elem1, elem2, max_ratio in constraints:
            print(f"  - {elem1}:{elem2} 比例 <= {max_ratio}:1")

        if exclude_element_counts:
            exclude_str = ", ".join(f"{n}元" for n in sorted(exclude_element_counts))
            print(f"  - 排除：{exclude_str}化合物不移动")

        print()

    # 统计所有文件
    all_compositions: Dict[Path, OrderedDict[str, int]] = {}
    total_files = 0

    for it_dir in it_dirs:
        res_files = sorted(it_dir.glob("*.res"))
        total_files += len(res_files)

    print(f"扫描目录：{', '.join(d.name for d in it_dirs)}")
    print(f"文件总数：{total_files}\n")

    # 收集所有文件路径
    all_files = []
    for it_dir in it_dirs:
        all_files.extend(sorted(it_dir.glob("*.res")))

    # 并行解析所有文件（提速）
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"使用 {num_workers} 个进程并行处理\n")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_path = {executor.submit(parse_composition, path): path for path in all_files}

        for idx, future in enumerate(as_completed(future_to_path), 1):
            res_path = future_to_path[future]
            try:
                composition = future.result()
                all_compositions[res_path] = composition
            except Exception as e:
                print(f"\n警告：解析 {res_path} 时出错: {e}")

            show_progress_bar(idx, total_files)

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 统计配比分布
    composition_counter: Counter[str] = Counter()
    valid_counter: Counter[str] = Counter()
    invalid_counter: Counter[str] = Counter()

    valid_count = 0
    invalid_count = 0

    for res_path, composition in all_compositions.items():
        signature = composition_signature(composition)
        composition_counter[signature] += 1

        if constraints:
            is_valid, _ = check_all_constraints(composition, constraints, exclude_element_counts)
            if is_valid:
                valid_count += 1
                valid_counter[signature] += 1
            else:
                invalid_count += 1
                invalid_counter[signature] += 1

    # 生成详细报告
    detailed_report = generate_detailed_report(
        it_dirs,
        total_files,
        all_compositions,
        composition_counter,
        constraints,
        valid_counter if constraints else None,
        invalid_counter if constraints else None
    )

    # 自动保存详细报告到文件（固定名称，覆盖旧文件）
    default_report_name = "composition_report.txt"
    report_path = Path(output_file) if output_file else Path(default_report_name)

    report_path.write_text(detailed_report, encoding="utf-8")
    print(f"详细报告已保存到：{report_path}")

    # 如果有筛选条件，额外生成不满足条件的文件列表（用于快速移动）
    invalid_list_path = None
    if constraints and invalid_count > 0:
        invalid_list_name = "invalid_structures.txt"
        invalid_list_path = Path(invalid_list_name)

        # 收集所有不满足条件的文件路径和配比信息（应用 exclude_element_counts 过滤）
        invalid_files_info: List[Tuple[Path, str, str]] = []
        for res_path, composition in all_compositions.items():
            is_valid, violations = check_all_constraints(composition, constraints, exclude_element_counts)
            if not is_valid:
                signature = composition_signature(composition)
                reason = "; ".join(violations)
                invalid_files_info.append((res_path, signature, reason))

        # 写入文件
        with invalid_list_path.open("w", encoding="utf-8") as f:
            f.write("# 不满足配比约束条件的结构文件列表\n")
            f.write(f"# 生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("#\n")
            f.write("# 筛选条件：\n")
            for elem1, elem2, max_ratio in constraints:
                f.write(f"#   - {elem1}:{elem2} 比例 <= {max_ratio}:1\n")
            if exclude_element_counts:
                exclude_str = ", ".join(f"{n}元" for n in sorted(exclude_element_counts))
                f.write(f"#   - 排除：{exclude_str}化合物不移动\n")
            f.write(f"#\n")
            f.write(f"# 不满足条件的文件数：{len(invalid_files_info)}\n")
            f.write("#\n")
            f.write("# 格式：文件路径 | 配比 | 违反原因\n")
            f.write("#" + "=" * 68 + "\n\n")

            for res_path, signature, reason in invalid_files_info:
                f.write(f"{res_path} | {signature} | {reason}\n")

        print(f"不满足条件的文件列表已保存到：{invalid_list_path}")

    print()

    # 输出简要统计到终端
    lines: List[str] = []
    lines.append("=" * 70)
    lines.append("统计结果（简要）")
    lines.append("=" * 70)
    lines.append(f"文件总数：{total_files}")

    if constraints:
        lines.append(f"满足条件：{valid_count} ({valid_count*100/total_files:.1f}%)")
        lines.append(f"不满足条件：{invalid_count} ({invalid_count*100/total_files:.1f}%)")
        lines.append("")

    # 显示元素范围摘要
    lines.append("")
    lines.append("元素数量范围：")
    element_stats = analyze_element_ranges(all_compositions)
    for elem in sorted(element_stats.keys()):
        stats = element_stats[elem]
        lines.append(f"  {elem}: {stats['min']} ~ {stats['max']} (平均: {stats['avg']:.2f})")

    lines.append("")
    lines.append("配比数量（前10个）：")
    for signature, count in composition_counter.most_common(10):
        percentage = count * 100 / total_files
        lines.append(f"  {signature:<20} {count:>5} ({percentage:>5.1f}%)")

    if len(composition_counter) > 10:
        lines.append(f"  ... (还有 {len(composition_counter) - 10} 种配比)")

    lines.append("")
    lines.append(f"完整报告请查看：{report_path}")
    lines.append("=" * 70)

    print("\n".join(lines))


def filter_mode(
    it_dirs: List[Path],
    constraints: List[Tuple[str, str, float]],
    output_file: Optional[str]
):
    """筛选模式：列出满足条件的文件（并行处理）"""
    if not constraints:
        print("错误：筛选模式需要指定 --ratio-limit 参数", file=sys.stderr)
        return

    print("\n" + "=" * 70)
    print("配比筛选模式")
    print("=" * 70 + "\n")

    print("筛选条件：")
    for elem1, elem2, max_ratio in constraints:
        print(f"  - {elem1}:{elem2} 比例 <= {max_ratio}:1")
    print()

    # 解析所有文件
    all_files = []
    for it_dir in it_dirs:
        all_files.extend(sorted(it_dir.glob("*.res")))

    total_files = len(all_files)
    print(f"扫描文件总数：{total_files}\n")

    # 并行解析所有文件
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"使用 {num_workers} 个进程并行处理\n")

    valid_files: List[Tuple[Path, str]] = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_path = {executor.submit(parse_composition, path): path for path in all_files}

        for idx, future in enumerate(as_completed(future_to_path), 1):
            res_path = future_to_path[future]
            try:
                composition = future.result()
                is_valid, _ = check_all_constraints(composition, constraints)

                if is_valid:
                    signature = composition_signature(composition)
                    valid_files.append((res_path, signature))
            except Exception as e:
                print(f"\n警告：解析 {res_path} 时出错: {e}")

            show_progress_bar(idx, total_files)

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 输出结果
    lines: List[str] = []
    lines.append("=" * 70)
    lines.append(f"满足条件的文件：{len(valid_files)} / {total_files}")
    lines.append("=" * 70)

    for res_path, signature in valid_files:
        lines.append(f"{res_path.parent.name}/{res_path.name:<40} {signature}")

    output_text = "\n".join(lines)

    if output_file:
        output_path = Path(output_file)
        output_path.write_text(output_text + "\n", encoding="utf-8")
        print(f"结果已写入：{output_path}")
    else:
        print(output_text)


def read_invalid_structures_from_file(file_path: Path) -> List[Tuple[Path, str, str]]:
    """
    从文件中读取不满足条件的结构列表

    返回：[(文件路径, 配比, 违反原因), ...]
    """
    if not file_path.exists():
        print(f"错误：文件不存在: {file_path}", file=sys.stderr)
        return []

    invalid_files: List[Tuple[Path, str, str]] = []

    try:
        with file_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                # 跳过注释和空行
                if not line or line.startswith("#"):
                    continue

                # 解析格式：文件路径 | 配比 | 违反原因
                parts = line.split("|")
                if len(parts) >= 3:
                    file_path_str = parts[0].strip()
                    signature = parts[1].strip()
                    reason = parts[2].strip()
                    invalid_files.append((Path(file_path_str), signature, reason))

        print(f"[INFO] 从 {file_path} 读取到 {len(invalid_files)} 个不满足条件的结构")
        return invalid_files

    except Exception as e:
        print(f"错误：读取文件 {file_path} 时出错: {e}", file=sys.stderr)
        return []


def move_invalid_mode(
    it_dirs: List[Path],
    constraints: List[Tuple[str, str, float]],
    trash_dir: Path,
    skip_confirm: bool,
    exclude_element_counts: Optional[List[int]] = None
):
    """移动不合格结构到垃圾箱（并行处理）"""
    if not constraints:
        print("错误：移动模式需要指定 --ratio-limit 参数", file=sys.stderr)
        return

    print("\n" + "=" * 70)
    print("配比筛选工具 - 移动不合格结构")
    print("=" * 70 + "\n")

    print("筛选条件：")
    for elem1, elem2, max_ratio in constraints:
        print(f"  - {elem1}:{elem2} 比例 <= {max_ratio}:1")

    if exclude_element_counts:
        exclude_str = ", ".join(f"{n}元" for n in sorted(exclude_element_counts))
        print(f"  - 排除：{exclude_str}化合物不移动")

    print()

    # 解析所有文件
    all_files = []
    for it_dir in it_dirs:
        all_files.extend(sorted(it_dir.glob("*.res")))

    total_files = len(all_files)
    print(f"扫描目录：{', '.join(d.name for d in it_dirs)}")
    print(f"文件总数：{total_files}\n")

    # 分类文件
    valid_files: List[Path] = []
    invalid_files: List[Tuple[Path, str, str]] = []  # (路径, 配比, 违反原因)

    for idx, res_path in enumerate(all_files, 1):
        composition = parse_composition(res_path)
        is_valid, violations = check_all_constraints(composition, constraints, exclude_element_counts)

        if is_valid:
            valid_files.append(res_path)
        else:
            signature = composition_signature(composition)
            reason = "; ".join(violations)
            invalid_files.append((res_path, signature, reason))

        show_progress_bar(idx, total_files)

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 显示统计
    print("=" * 70)
    print("统计结果：")
    print(f"  文件总数：{total_files}")
    print(f"  满足条件：{len(valid_files)} ({len(valid_files)*100/total_files:.1f}%)")
    print(f"  不满足条件：{len(invalid_files)} ({len(invalid_files)*100/total_files:.1f}%)")
    print("=" * 70 + "\n")

    if not invalid_files:
        print("所有文件都满足条件，无需移动。")
        return

    # 显示不满足条件的配比分布
    invalid_signatures: Counter[str] = Counter()
    for _, signature, _ in invalid_files:
        invalid_signatures[signature] += 1

    print("不满足条件的配比分布（前10个）：")
    for signature, count in invalid_signatures.most_common(10):
        print(f"  {signature:<20} {count}")
    print()

    # 显示将要移动的文件（前20个）
    print(f"将要移动的文件（显示前20个，共{len(invalid_files)}个）：")
    for res_path, signature, reason in invalid_files[:20]:
        it_name = res_path.parent.name
        print(f"  {it_name}/{res_path.name:<35} {signature:<15} ({reason})")

    if len(invalid_files) > 20:
        print(f"  ... (还有 {len(invalid_files) - 20} 个文件)")
    print()

    # 确认操作
    print(f"⚠️  警告：即将移动 {len(invalid_files)} 个文件到 {trash_dir}/")

    if not skip_confirm:
        try:
            response = input("确认继续？[y/N]: ").strip().lower()
            if response not in ['y', 'yes']:
                print("操作已取消。")
                return
        except (KeyboardInterrupt, EOFError):
            print("\n操作已取消。")
            return
    print()

    # 并行执行移动
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    moved_files: List[Tuple[Path, Path, str]] = []
    failed_files: List[Tuple[Path, str]] = []

    # 准备移动任务
    move_tasks = []
    for res_path, signature, reason in invalid_files:
        it_name = res_path.parent.name
        move_tasks.append((res_path, trash_dir, it_name, reason))

    # 使用多进程并行移动
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"使用 {num_workers} 个进程并行移动文件\n")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(move_single_file, task): task for task in move_tasks}

        for idx, future in enumerate(as_completed(futures), 1):
            try:
                success, src_path, target_or_error, reason = future.result()

                if success:
                    moved_files.append((src_path, target_or_error, reason))
                else:
                    failed_files.append((src_path, target_or_error))
            except Exception as e:
                task = futures[future]
                failed_files.append((task[0], str(e)))

            if sys.stdout.isatty():
                filled = int(40 * idx / len(invalid_files))
                bar = "#" * filled + "." * (40 - filled)
                sys.stdout.write(f"\r[移动] [{bar}] {idx}/{len(invalid_files)}")
                sys.stdout.flush()

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 写入日志（固定名称，覆盖旧文件）
    if moved_files:
        log_path = trash_dir / "move_log.txt"
        write_move_log(log_path, moved_files, constraints)

    # 显示结果
    print("=" * 70)
    print("操作完成！")
    print(f"  成功移动：{len(moved_files)} 个文件")
    if moved_files:
        print(f"  操作日志：{trash_dir / 'move_log.txt'}")

    if failed_files:
        print(f"  移动失败：{len(failed_files)} 个文件")
        for res_path, error in failed_files[:5]:
            print(f"    {res_path}: {error}")
        if len(failed_files) > 5:
            print(f"    ... (还有 {len(failed_files) - 5} 个)")

    print("\nPD 目录状态：")
    print(f"  剩余有效结构：{len(valid_files)} 个")
    print(f"  已移动到垃圾箱：{len(moved_files)} 个")
    print("=" * 70)


def restore_from_trash(
    trash_dir: Path,
    log_file: Optional[Path],
    skip_confirm: bool
):
    """一键还原：从垃圾箱还原文件到原位置（并行处理）"""
    print("\n" + "=" * 70)
    print("配比筛选工具 - 一键还原")
    print("=" * 70 + "\n")

    if not trash_dir.exists():
        print(f"错误：垃圾箱目录不存在: {trash_dir}")
        return

    # 如果指定了日志文件，从日志读取
    if log_file and log_file.exists():
        print(f"从日志文件读取移动记录：{log_file}\n")

        restore_tasks: List[Tuple[Path, Path]] = []

        try:
            with log_file.open("r", encoding="utf-8") as f:
                in_details = False
                for line in f:
                    line = line.strip()
                    if line.startswith("移动详情"):
                        in_details = True
                        continue
                    if not in_details:
                        continue
                    if line.startswith("源文件："):
                        src_path_str = line.replace("源文件：", "").strip()
                    elif line.startswith("目标："):
                        target_path_str = line.replace("目标：", "").strip()
                        # 创建还原任务：从target还原到src
                        src_path = Path(src_path_str)
                        target_path = Path(target_path_str)
                        if target_path.exists():
                            restore_tasks.append((target_path, src_path))
        except Exception as e:
            print(f"警告：读取日志文件失败: {e}")
            print("将扫描整个垃圾箱目录...\n")
            restore_tasks = []

    # 如果没有日志或读取失败，扫描垃圾箱
    if not restore_tasks:
        print("扫描垃圾箱目录...\n")

        for it_dir in trash_dir.glob("IT*"):
            if not it_dir.is_dir():
                continue

            it_name = it_dir.name
            for trash_file in it_dir.glob("*.res"):
                # 尝试还原到原位置
                original_path = Path(it_name) / trash_file.name
                restore_tasks.append((trash_file, original_path))

    if not restore_tasks:
        print("垃圾箱为空，没有需要还原的文件。")
        return

    print(f"找到 {len(restore_tasks)} 个文件需要还原\n")

    # 显示前20个文件
    print("将要还原的文件（显示前20个）：")
    for trash_path, original_path in restore_tasks[:20]:
        print(f"  {trash_path} → {original_path}")

    if len(restore_tasks) > 20:
        print(f"  ... (还有 {len(restore_tasks) - 20} 个文件)")
    print()

    # 确认操作
    print(f"⚠️  警告：即将还原 {len(restore_tasks)} 个文件")

    if not skip_confirm:
        try:
            response = input("确认继续？[y/N]: ").strip().lower()
            if response not in ['y', 'yes']:
                print("操作已取消。")
                return
        except (KeyboardInterrupt, EOFError):
            print("\n操作已取消。")
            return
    print()

    # 并行执行还原
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    restored_files: List[Tuple[Path, Path]] = []
    failed_files: List[Tuple[Path, str]] = []

    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"使用 {num_workers} 个进程并行还原文件\n")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(restore_single_file, task): task for task in restore_tasks}

        for idx, future in enumerate(as_completed(futures), 1):
            try:
                success, trash_path, target_or_error, _ = future.result()

                if success:
                    restored_files.append((trash_path, target_or_error))
                else:
                    failed_files.append((trash_path, target_or_error))
            except Exception as e:
                task = futures[future]
                failed_files.append((task[0], str(e)))

            if sys.stdout.isatty():
                filled = int(40 * idx / len(restore_tasks))
                bar = "#" * filled + "." * (40 - filled)
                sys.stdout.write(f"\r[还原] [{bar}] {idx}/{len(restore_tasks)}")
                sys.stdout.flush()

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 显示结果
    print("=" * 70)
    print("还原完成！")
    print(f"  成功还原：{len(restored_files)} 个文件")

    if failed_files:
        print(f"  还原失败：{len(failed_files)} 个文件")
        for trash_path, error in failed_files[:5]:
            print(f"    {trash_path}: {error}")
        if len(failed_files) > 5:
            print(f"    ... (还有 {len(failed_files) - 5} 个)")

    print("=" * 70)


def quick_move_mode(
    invalid_list_file: Path,
    trash_dir: Path,
    skip_confirm: bool,
    dry_run: bool = False
):
    """快速移动模式：直接从文件列表移动，不重新解析（并行处理）"""
    print("\n" + "=" * 70)
    print("配比筛选工具 - 快速移动模式（读取文件列表）")
    print("=" * 70 + "\n")

    # 读取文件列表
    invalid_files = read_invalid_structures_from_file(invalid_list_file)

    if not invalid_files:
        print("错误：无法读取到任何文件")
        return

    total_count = len(invalid_files)

    # 检查文件存在性
    existing_files: List[Tuple[Path, str, str]] = []
    missing_files: List[Path] = []

    for res_path, signature, reason in invalid_files:
        if res_path.exists():
            existing_files.append((res_path, signature, reason))
        else:
            missing_files.append(res_path)

    print(f"文件总数：{total_count}")
    print(f"存在的文件：{len(existing_files)}")
    if missing_files:
        print(f"不存在的文件：{len(missing_files)}")
        print()

    if not existing_files:
        print("错误：所有文件都不存在，无法移动")
        return

    # 显示配比分布
    signature_counter: Counter[str] = Counter()
    for _, signature, _ in existing_files:
        signature_counter[signature] += 1

    print("不满足条件的配比分布（前10个）：")
    for signature, count in signature_counter.most_common(10):
        print(f"  {signature:<20} {count}")
    print()

    # 显示将要移动的文件（前20个）
    print(f"将要移动的文件（显示前20个，共{len(existing_files)}个）：")
    for res_path, signature, reason in existing_files[:20]:
        it_name = res_path.parent.name if res_path.parent.name.startswith("IT") else "unknown"
        print(f"  {it_name}/{res_path.name:<35} {signature:<15} ({reason})")

    if len(existing_files) > 20:
        print(f"  ... (还有 {len(existing_files) - 20} 个文件)")
    print()

    # 确认操作
    if dry_run:
        print("[模拟模式] 不会实际移动文件")
        print()
        return
    else:
        print(f"⚠️  警告：即将移动 {len(existing_files)} 个文件到 {trash_dir}/")

        if not skip_confirm:
            try:
                response = input("确认继续？[y/N]: ").strip().lower()
                if response not in ['y', 'yes']:
                    print("操作已取消。")
                    return
            except (KeyboardInterrupt, EOFError):
                print("\n操作已取消。")
                return
        print()

    # 并行执行移动
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, as_completed

    moved_files: List[Tuple[Path, Path, str]] = []
    failed_files: List[Tuple[Path, str]] = []

    # 准备移动任务
    move_tasks = []
    for res_path, signature, reason in existing_files:
        it_name = res_path.parent.name if res_path.parent.name.startswith("IT") else "unknown"
        move_tasks.append((res_path, trash_dir, it_name, reason))

    # 使用多进程并行移动
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"使用 {num_workers} 个进程并行移动文件\n")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(move_single_file, task): task for task in move_tasks}

        for idx, future in enumerate(as_completed(futures), 1):
            try:
                success, src_path, target_or_error, reason = future.result()

                if success:
                    moved_files.append((src_path, target_or_error, reason))
                else:
                    failed_files.append((src_path, target_or_error))
            except Exception as e:
                task = futures[future]
                failed_files.append((task[0], str(e)))

            if sys.stdout.isatty():
                filled = int(40 * idx / len(existing_files))
                bar = "#" * filled + "." * (40 - filled)
                sys.stdout.write(f"\r[移动] [{bar}] {idx}/{len(existing_files)}")
                sys.stdout.flush()

    if sys.stdout.isatty():
        sys.stdout.write("\n\n")

    # 写入日志（固定名称，覆盖旧文件）
    if moved_files:
        log_path = trash_dir / "move_log.txt"
        write_move_log(log_path, moved_files, [])  # 空约束列表

    # 显示结果
    print("=" * 70)
    print("操作完成！")
    print(f"  成功移动：{len(moved_files)} 个文件")
    if moved_files:
        print(f"  操作日志：{trash_dir / 'move_log.txt'}")

    if failed_files:
        print(f"  移动失败：{len(failed_files)} 个文件")
        for res_path, error in failed_files[:5]:
            print(f"    {res_path}: {error}")
        if len(failed_files) > 5:
            print(f"    ... (还有 {len(failed_files) - 5} 个)")

    if missing_files:
        print(f"  文件不存在（跳过）：{len(missing_files)} 个")

    print("=" * 70)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="配比筛选工具：统计、筛选、移动不满足配比约束条件的结构",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 统计单个IT目录的配比分布
  python filter_composition.py IT0 --stats

  # 统计所有IT*目录，并应用筛选条件
  python filter_composition.py --ratio-limit Ce:Mg:5 --stats

  # 筛选满足条件的文件
  python filter_composition.py IT0 --ratio-limit Ce:Mg:5 --filter

  # 移动所有不满足条件的文件到垃圾箱（交互确认）
  python filter_composition.py --ratio-limit Ce:Mg:5 --move-invalid

  # 多个约束条件
  python filter_composition.py --ratio-limit Ce:Mg:5 --ratio-limit H:Mg:3 --move-invalid

  # 跳过确认直接移动
  python filter_composition.py --ratio-limit Ce:Mg:5 --move-invalid --yes

  # 排除单质、二元、三元化合物，只移动四元及以上化合物
  python filter_composition.py --ratio-limit Ce:Mg:5 --move-invalid --exclude-elements 1 2 3

  # 快速移动模式（推荐工作流）
  # 步骤1：先统计，生成报告和文件列表
  python filter_composition.py --ratio-limit Ce:Mg:5 --stats
  # 生成: composition_report_*.txt 和 invalid_structures_*.txt

  # 步骤2：查看报告后，直接用文件列表快速移动（不重新解析）
  python filter_composition.py --read-file invalid_structures.txt

  # 快速移动 + 模拟运行（仅快速移动模式支持 --dry-run）
  python filter_composition.py --read-file invalid_structures.txt --dry-run

  # 一键还原：将垃圾箱中的所有文件还原
  python filter_composition.py --restore

  # 从指定日志文件还原
  python filter_composition.py --restore --restore-log trash_bin/move_log.txt

  # 跳过确认直接还原
  python filter_composition.py --restore --yes

注意事项:
  - 需要在 PD 目录内运行
  - 默认扫描所有 IT* 目录，可通过 ITERATION 参数指定单个目录
  - 比例约束格式：ELEM1:ELEM2:MAX_RATIO（例如 Ce:Mg:5 表示 Ce:Mg <= 5:1 或 1:5）
  - 移动的文件会保留 IT* 目录结构
  - --stats 模式会自动生成不满足条件的文件列表，可用于快速移动
  - 标准移动模式不支持 --dry-run，仅快速移动模式支持
  - --exclude-elements 用于排除特定元素数量的化合物（如单质、二元化合物等）
        """
    )

    parser.add_argument(
        "iteration",
        nargs="?",
        help="迭代编号（可选），例如 0 或 IT0。不指定则扫描所有 IT* 目录。",
    )

    parser.add_argument(
        "--ratio-limit",
        action="append",
        metavar="ELEM1:ELEM2:RATIO",
        help="元素比例约束，格式：ELEM1:ELEM2:MAX_RATIO。可多次指定。"
             "例如 Ce:Mg:5 表示 Ce 和 Mg 的比例不超过 1:5 或 5:1。",
    )

    # 模式选择
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--stats",
        action="store_true",
        help="统计模式（默认）：显示配比分布统计",
    )
    mode_group.add_argument(
        "--filter",
        action="store_true",
        help="筛选模式：列出满足条件的文件",
    )
    mode_group.add_argument(
        "--move-invalid",
        action="store_true",
        help="移动模式：将不满足条件的文件移动到垃圾箱",
    )
    mode_group.add_argument(
        "--restore",
        action="store_true",
        help="一键还原模式：将trash_bin中的文件还原到原位置",
    )

    # 移动选项
    parser.add_argument(
        "--trash-dir",
        default="trash_bin",
        help="垃圾箱目录（默认：PD/trash_bin）",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="[仅用于快速移动模式] 模拟运行：显示将要移动的文件，但不实际移动",
    )
    parser.add_argument(
        "--yes", "-y",
        action="store_true",
        help="跳过确认，直接执行移动（危险操作）",
    )
    parser.add_argument(
        "--read-file",
        metavar="FILE",
        help="快速移动模式：直接从文件列表读取并移动（跳过重新解析所有结构）。"
             "通常使用 --stats 生成的 invalid_structures_*.txt 文件。支持 --dry-run 选项。",
    )
    parser.add_argument(
        "--restore-log",
        metavar="FILE",
        help="[用于还原模式] 指定还原时使用的日志文件（可选，不指定则扫描整个垃圾箱）",
    )
    parser.add_argument(
        '--exclude-elements',
        type=int,
        nargs='+',
        metavar='N',
        help='指定不移动的化合物元素数量。例如：1 2 3 表示单质、二元、三元化合物不移动',
    )

    # 输出选项
    parser.add_argument(
        "--output",
        help="保存结果到文件",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    # 检查当前目录
    pd_dir = Path.cwd()
    if pd_dir.name != "PD":
        print("警告：当前不在 PD 目录，建议切换到 PD 再运行。", file=sys.stderr)

    # 解析迭代目录
    if args.iteration:
        try:
            iteration = parse_iteration_label(args.iteration)
        except ValueError as exc:
            print(f"错误：{exc}", file=sys.stderr)
            return 1
        it_dirs = find_it_directories(pd_dir, iteration)
        if not it_dirs:
            print(f"错误：未找到迭代目录 {iteration}", file=sys.stderr)
            return 1
    else:
        it_dirs = find_it_directories(pd_dir)
        if not it_dirs:
            print("错误：未找到任何 IT* 目录", file=sys.stderr)
            return 1

    # 解析比例约束
    constraints: List[Tuple[str, str, float]] = []
    if args.ratio_limit:
        for ratio_str in args.ratio_limit:
            try:
                constraint = parse_ratio_limit(ratio_str)
                constraints.append(constraint)
            except ValueError as exc:
                print(f"错误：{exc}", file=sys.stderr)
                return 1

    # 确定运行模式
    if args.restore:
        # 还原模式
        print("[INFO] 运行模式: 一键还原模式")
        print()

        trash_dir = pd_dir / args.trash_dir
        log_file = Path(args.restore_log) if args.restore_log else None
        restore_from_trash(trash_dir, log_file, args.yes)
    elif args.read_file:
        # 快速移动模式：直接从文件读取
        print("[INFO] 运行模式: 快速移动模式（读取文件列表）")
        print()

        invalid_list_file = Path(args.read_file)
        trash_dir = pd_dir / args.trash_dir
        quick_move_mode(invalid_list_file, trash_dir, args.yes, args.dry_run)
    elif args.move_invalid:
        # 标准移动模式：重新解析所有结构
        print("[INFO] 运行模式: 标准移动模式（重新解析所有结构）")

        # 检查是否使用了 --dry-run（标准移动模式不支持）
        if args.dry_run:
            print("错误：标准移动模式不支持 --dry-run 选项。", file=sys.stderr)
            print("如需模拟运行，请使用快速移动模式：", file=sys.stderr)
            print("  1. 先运行 --stats 生成文件列表", file=sys.stderr)
            print("  2. 然后使用 --read-file <文件列表> --dry-run", file=sys.stderr)
            return 1

        print()

        trash_dir = pd_dir / args.trash_dir
        exclude_counts = args.exclude_elements if args.exclude_elements else None
        move_invalid_mode(it_dirs, constraints, trash_dir, args.yes, exclude_counts)
    elif args.filter:
        # 筛选模式
        filter_mode(it_dirs, constraints, args.output)
    else:
        # 默认统计模式
        exclude_counts = args.exclude_elements if args.exclude_elements else None
        stats_mode(it_dirs, constraints, args.output, exclude_counts)

    return 0


if __name__ == "__main__":
    sys.exit(main())
