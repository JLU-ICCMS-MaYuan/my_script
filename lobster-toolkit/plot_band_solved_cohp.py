# -*- coding: utf-8 -*-
"""
这个脚本用于解析LOBSTER输出的一系列KspaceCOHPBand*.lobster文件,
并将它们合并，用于绘制能带分辨的投影晶体轨道哈密顿布居 (band-resolved projected COHP) 图。

核心功能特点:
- 自动查找并按序号排序处理所有 KspaceCOHPBand*.lobster 文件
- 支持多种轨道相互作用的同时分析（通过正则表达式匹配）
- 读取POSCAR文件计算倒易晶格，精确计算K点路径距离
- 从KPOINTS文件读取高对称路径上的K点分数坐标
- 从KPATH.in文件获取高对称点标签并在图中标记
- 输出.dat文件包含：K点距离、能量、各种pCOHP值（多列格式）
- 输出.png图片：胖能带风格的散点图，颜色表示总pCOHP贡献

输出文件说明:
- .dat文件: 第一列为K点距离，第二列为能量，后续列为各pCOHP值
- .dat文件按能带分块，每个能带的数据后有空行分隔
- .png图片: 背景为灰色能带线，散点颜色表示-pCOHP值的总和

如何使用:
1. 准备 VASP 计算生成的 POSCAR, KPOINTS 和 KPATH.in 文件
2. 确保所有 KspaceCOHPBand*.lobster 文件在同一目录下
3. 修改下方"用户设置"部分的 `interactions_to_plot` 列表，定义要分析的轨道相互作用
4. 根据需要调整费米能级和绘图参数
5. 运行脚本: python plot_kcohp.py

注意事项:
- 确保轨道标签的正则表达式与LOBSTER输出文件中的标签完全匹配
- 可以通过查看KspaceCOHPBand1.lobster文件的第3行来确认可用的轨道标签
"""

import re
import glob
import numpy as np
import matplotlib.pyplot as plt

# --- 用户设置 ---

# 1. 定义您想要分析的一个或多个相互作用 (化学键)。
#    这是一个列表，您可以放入任意数量的字典，每个字典代表一种相互作用。
#    - 'label': 图例、文件列名以及输出文件名中显示的标签。
#    - 'regex1', 'regex2': 可以是单个字符串(匹配一类轨道)，也可以是字符串列表(精确匹配多个轨道)。
#    注意：根据LOBSTER输出文件中的实际轨道标签进行匹配
interactions_to_plot = [
    {
        "label": "La_4f_A2u", 
        "regex1": [r"La1_4f_z\^3"],
        "regex2": r"H\d+_1s"
    },
    {
        "label": "La_4f_B2u",
        "regex1": [r"La1_4f_x\(x\^2-3y\^2\)"],
        "regex2": r"H\d+_1s"
    },
    {
        "label": "La_4f_E1u",
        "regex1": [r"La1_4f_xz\^2", r"La1_4f_yz\^2"],
        "regex2": r"H\d+_1s"
    },
    {
        "label": "La_4f_E2u",
        "regex1": [r"La1_4f_xyz", r"La1_4f_z\(x\^2-y\^2\)"],
        "regex2": r"H\d+_1s"
    },
    {
        "label": "La_4f_B1u",
        "regex1": [r"La1_4f_y\(3x\^2-y\^2\)"],
        "regex2": r"H\d+_1s"
    },
    {
        "label": "La_4f_all_H",
        "regex1": r"La\d+_4f_.*",
        "regex2": r"H\d+_1s"
    },
    # 如果想分析Sc的轨道相互作用，可以添加如下示例：
    # {
    #     "label": "Sc_3d_H",
    #     "regex1": r"Sc\d+_3d_.*",
    #     "regex2": r"H\d+_1s"
    # },
]


# 2. 绘图参数设置
FERMI_ENERGY = 12.5967  # 费米能级，能量坐标会相对于此值进行偏移。请根据lobsterout文件确认费米能级
ENERGY_RANGE = [-10, 5]  # 绘图的能量范围 [E_min, E_max] (相对于费米能级)
VMAX = 1.0  # 颜色映射的最大值 (-pCOHP)。可以调整这个值来改变颜色对比度

# 提示：如何确定合适的参数值
# - FERMI_ENERGY: 查看lobsterout文件中的"Fermi energy"行
# - ENERGY_RANGE: 根据感兴趣的能量窗口调整，通常价带和导带边缘附近最重要
# - VMAX: 根据实际pCOHP值的范围调整，可以先运行一次查看数据范围再优化

# --- 核心功能代码 ---

def get_reciprocal_lattice_from_poscar(poscar_file="POSCAR"):
    """
    从POSCAR文件读取晶格矢量并计算倒易晶格矢量。
    """
    try:
        with open(poscar_file, 'r') as f:
            lines = f.readlines()
        
        scale = float(lines[1].strip())
        a1 = np.array(list(map(float, lines[2].split()))) * scale
        a2 = np.array(list(map(float, lines[3].split()))) * scale
        a3 = np.array(list(map(float, lines[4].split()))) * scale
        
        volume = np.dot(a1, np.cross(a2, a3))
        
        b1 = 2 * np.pi * np.cross(a2, a3) / volume
        b2 = 2 * np.pi * np.cross(a3, a1) / volume
        b3 = 2 * np.pi * np.cross(a1, a2) / volume
        
        return np.array([b1, b2, b3])

    except FileNotFoundError:
        print(f"错误: 未找到 {poscar_file} 文件。无法计算倒易晶格，程序终止。")
        return None
    except (IOError, IndexError, ValueError) as e:
        print(f"错误: 解析 {poscar_file} 文件时出错: {e}")
        return None

def parse_kpoints_and_kpath(kpoints_file="KPOINTS", kpath_file="KPATH.in", poscar_file="POSCAR"):
    """
    解析 KPOINTS, KPATH.in 和 POSCAR 文件。
    - 从 POSCAR 计算倒易晶格。
    - 从 KPATH.in 获取高对称点的坐标和标签。
    - 从 KPOINTS 获取高对称路径上的K点分数坐标。
    - 计算沿路径的笛卡尔坐标累积距离，同时识别高对称点。
    """
    reciprocal_lattice = get_reciprocal_lattice_from_poscar(poscar_file)
    if reciprocal_lattice is None:
        return None, None

    # 先解析 KPATH.in 获取高对称点信息
    high_sym_points_dict = {}
    try:
        with open(kpath_file, 'r') as f:
            lines = f.readlines()
        
        # 过滤掉标题行和注释行，只保留包含坐标和标签的行
        for line in lines:
            line = line.strip()
            if not line or "!" in line or "K-Path" in line or "Line-Mode" in line or "Reciprocal" in line:
                continue
            parts = line.split()
            # 检查是否是坐标行（至少4列：x, y, z, label）
            if len(parts) >= 4:
                try:
                    # 测试前三列是否为数字
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                    label = parts[-1]
                    if "GAMMA" in label.upper():
                        label = r"$\Gamma$"
                    # 使用坐标作为键，避免重复
                    coord_key = (x, y, z)
                    high_sym_points_dict[coord_key] = label
                except ValueError:
                    continue
        
        print(f"成功从 {kpath_file} 读取 {len(high_sym_points_dict)} 个高对称点定义。")

    except FileNotFoundError:
        print(f"警告: 未找到 {kpath_file} 文件。将不标记高对称点。")
        high_sym_points_dict = {}
    except Exception as e:
        print(f"警告: 解析 {kpath_file} 时出错: {e}。将不标记高对称点。")
        high_sym_points_dict = {}

    # 解析 KPOINTS 文件
    try:
        with open(kpoints_file, 'r') as f:
            lines = f.readlines()
        
        parts = lines[0].split()
        if len(parts) < 7:
            print(f"错误: {kpoints_file} 文件的第一行格式不正确。")
            return None, None
            
        num_ibz_kpts = int(parts[4])
        num_path_kpts = int(parts[6])

        path_kpoints_coords = []
        start_line_index = 3 + num_ibz_kpts
        
        end_line_index = start_line_index + num_path_kpts
        if end_line_index > len(lines):
            end_line_index = len(lines)

        for i in range(start_line_index, end_line_index):
            try:
                coords = list(map(float, lines[i].split()[:3]))
                path_kpoints_coords.append(coords)
            except (ValueError, IndexError):
                continue
        
        actual_path_kpts = len(path_kpoints_coords)
        if actual_path_kpts != num_path_kpts:
             print(f"警告: 最终成功解析 {actual_path_kpts} 个路径K点，与文件头声明的 {num_path_kpts} 个不符。")
             num_path_kpts = actual_path_kpts

        print(f"成功从 {kpoints_file} 读取 {num_path_kpts} 个高对称路径K点。")

    except FileNotFoundError:
        print(f"错误: 未找到 {kpoints_file} 文件。程序终止。")
        return None, None
    except Exception as e:
        print(f"错误: 解析 {kpoints_file} 文件时出错: {e}")
        return None, None

    # 计算累积距离 (笛卡尔坐标) 并同时识别高对称点
    distances = [0.0]
    high_sym_points = []
    
    # 检查第一个点是否为高对称点
    if path_kpoints_coords:
        first_coord = tuple(path_kpoints_coords[0])
        for sym_coord, sym_label in high_sym_points_dict.items():
            if np.allclose(first_coord, sym_coord, atol=1e-4):
                high_sym_points.append({'label': sym_label, 'distance': 0.0})
                break
    
    # 计算累积距离并识别高对称点
    for i in range(1, len(path_kpoints_coords)):
        prev_k_frac = np.array(path_kpoints_coords[i-1])
        curr_k_frac = np.array(path_kpoints_coords[i])
        prev_k_cart = np.dot(prev_k_frac, reciprocal_lattice)
        curr_k_cart = np.dot(curr_k_frac, reciprocal_lattice)
        dist = np.linalg.norm(curr_k_cart - prev_k_cart)
        current_distance = distances[-1] + dist
        distances.append(current_distance)
        
        # 检查当前点是否为高对称点
        curr_coord = tuple(path_kpoints_coords[i])
        for sym_coord, sym_label in high_sym_points_dict.items():
            if np.allclose(curr_coord, sym_coord, atol=1e-4):
                high_sym_points.append({'label': sym_label, 'distance': current_distance})
                break
    
    # 去除重复的高对称点（基于距离）
    unique_high_sym_points = []
    seen_distances = set()
    for point in sorted(high_sym_points, key=lambda x: x['distance']):
        # 使用较小的容差来避免重复
        if not any(abs(point['distance'] - seen_dist) < 1e-6 for seen_dist in seen_distances):
            unique_high_sym_points.append(point)
            seen_distances.add(point['distance'])

    print(f"成功识别 {len(unique_high_sym_points)} 个高对称点在路径上。")
    return np.array(distances), unique_high_sym_points


def parse_single_band_cohp_file(filename):
    """解析单个 KspaceCOHPBand<N>.lobster 文件。"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError: return None
    try:
        header_match = re.search(r'kpoints\s+(\d+)', lines[0])
        num_kpoints = int(header_match.group(1))
    except (IndexError, ValueError, AttributeError): return None
    data_blocks = "".join(lines).split('COHP_munujk matrix for spin')[1:]
    if len(data_blocks) != num_kpoints: num_kpoints = len(data_blocks)
    if not data_blocks: return None
    labels = data_blocks[0].split('\n')[1].split()[1:]
    num_orbitals = len(labels)
    energies_one_band = np.zeros(num_kpoints)
    cohps_one_band = np.zeros((num_kpoints, num_orbitals, num_orbitals))
    for ik, block in enumerate(data_blocks):
        block_lines = block.split('\n')
        eigenval_match = re.search(r'eigenvalue\s+([-\d.]+)', block_lines[0])
        if eigenval_match: energies_one_band[ik] = float(eigenval_match.group(1))
        matrix_data = [line.split()[1:] for line in block_lines[2:2 + num_orbitals] if line]
        if len(matrix_data) == num_orbitals:
            cohps_one_band[ik, :, :] = np.array(matrix_data, dtype=float)
    return num_kpoints, energies_one_band, cohps_one_band, labels

def calculate_pcohp(cohps, labels, regex1, regex2):
    """根据一个或多个正则表达式计算投影COHP (pCOHP)。"""
    def get_indices_from_regex(regex_input):
        regex_list = regex_input if isinstance(regex_input, (list, tuple)) else [regex_input]
        indices = set()
        for r in regex_list:
            indices.update([i for i, label in enumerate(labels) if re.fullmatch(r, label)])
        return sorted(list(indices))

    indices1 = get_indices_from_regex(regex1)
    indices2 = get_indices_from_regex(regex2)
    
    if not indices1 or not indices2:
        if not hasattr(calculate_pcohp, "warned_pairs"): calculate_pcohp.warned_pairs = set()
        pair_key = (str(regex1), str(regex2))
        if pair_key not in calculate_pcohp.warned_pairs:
            print(f"警告: 对于相互作用 '{str(regex1)}' 和 '{str(regex2)}'，未找到完全匹配的轨道对。")
            calculate_pcohp.warned_pairs.add(pair_key)
        return np.zeros(cohps.shape[0])

    pcohp_values = np.zeros(cohps.shape[0])
    for i in indices1:
        for j in indices2:
            pcohp_values += cohps[:, i, j]
    return pcohp_values

def save_data_to_dat(filename, k_distances, energies, pcohps_dict, fermi_energy=0.0):
    """将能带和多个pCOHP数据保存为.dat文件, 不同能带间用空行分隔。"""
    num_bands, num_kpoints = energies.shape
    shifted_energies = energies - fermi_energy
    labels = list(pcohps_dict.keys())
    
    with open(filename, 'w') as f:
        header_labels = ["# K-Distance", "Energy"] + labels
        f.write("   ".join(f"{lbl:<15}" for lbl in header_labels) + "\n")
        
        for ib in range(num_bands):
            for ik in range(num_kpoints):
                pcohp_vals_str = " ".join([f"{pcohps_dict[lbl][ib, ik]:<15.6f}" for lbl in labels])
                f.write(f"  {k_distances[ik]:<15.6f} {shifted_energies[ib, ik]:<15.6f} {pcohp_vals_str}\n")
            f.write("\n")
            
    print(f"\n数据已保存至: {filename}")

def calculate_subplot_layout(n_plots):
    """根据子图数量计算最佳的行列布局"""
    if n_plots == 1:
        return 1, 1
    elif n_plots == 2:
        return 1, 2
    elif n_plots <= 4:
        return 2, 2
    elif n_plots <= 6:
        return 2, 3
    elif n_plots <= 9:
        return 3, 3
    elif n_plots <= 12:
        return 3, 4
    elif n_plots <= 16:
        return 4, 4
    else:
        # 对于更多子图，计算接近正方形的布局
        import math
        cols = math.ceil(math.sqrt(n_plots))
        rows = math.ceil(n_plots / cols)
        return rows, cols

def plot_band_resolved_pcohp(filename, k_distances, energies, pcohps_dict, high_sym_points, fermi_energy=0.0, energy_range=None, vmax=1.0):
    """绘制 band-resolved pCOHP 多子图，每个相互作用单独显示"""
    num_bands, num_kpoints = energies.shape
    shifted_energies = energies - fermi_energy
    
    # 获取相互作用列表
    interaction_labels = list(pcohps_dict.keys())
    n_interactions = len(interaction_labels)
    
    # 计算子图布局
    rows, cols = calculate_subplot_layout(n_interactions)
    
    # 创建图形，调整尺寸以适应多子图
    fig_width = min(20, cols * 6)  # 每个子图宽度约6英寸，最大20英寸
    fig_height = min(16, rows * 4)  # 每个子图高度约4英寸，最大16英寸
    fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height))
    
    # 确保axes是二维数组，即使只有一个子图
    if n_interactions == 1:
        axes = np.array([[axes]])
    elif rows == 1 or cols == 1:
        axes = axes.reshape(rows, cols)
    
    # 设置全局颜色范围，确保所有子图使用相同的colorbar范围
    global_vmin = -vmax
    global_vmax = vmax
    
    # 为每个相互作用创建子图
    scatter_plots = []
    for i, label in enumerate(interaction_labels):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        
        # 绘制背景能带结构（灰色线条）
        for ib in range(num_bands):
            ax.plot(k_distances, shifted_energies[ib, :], color='grey', lw=0.5, alpha=0.7)
        
        # 获取当前相互作用的pCOHP数据
        current_pcohp = pcohps_dict[label]
        
        # 准备散点图数据
        x_coords = np.tile(k_distances, num_bands)
        y_coords = shifted_energies.flatten()
        colors = current_pcohp.flatten()
        
        # 绘制散点图
        sc = ax.scatter(x_coords, y_coords, c=colors, cmap='coolwarm_r', s=15, 
                       vmin=global_vmin, vmax=global_vmax, rasterized=True)
        scatter_plots.append(sc)
        
        # 设置坐标轴范围
        ax.set_xlim(k_distances[0], k_distances[-1])
        if energy_range: 
            ax.set_ylim(energy_range[0], energy_range[1])
        
        # 添加费米能级线
        ax.axhline(0, color='k', linestyle='--', linewidth=1, alpha=0.8)
        
        # 设置标题
        ax.set_title(f"{label}", fontsize=12, fontweight='bold')
        
        # 设置高对称点标记
        if high_sym_points:
            tick_locs = [p['distance'] for p in high_sym_points]
            tick_lbls = [p['label'] for p in high_sym_points]
            ax.set_xticks(tick_locs)
            # 只在最底行显示标签
            if row == rows - 1:
                ax.set_xticklabels(tick_lbls, fontsize=10)
            else:
                ax.set_xticklabels([])
            # 绘制高对称点的垂直线
            for loc in tick_locs:
                ax.axvline(x=loc, color='k', linestyle='--', linewidth=0.6, alpha=0.6)
        else:
            if row == rows - 1:
                ax.set_xlabel("K-point Path", fontsize=12)
            ax.set_xticks([])
        
        # 只在最左列显示y轴标签
        if col == 0:
            ax.set_ylabel(r"Energy (E - E$_F$) [eV]", fontsize=12)
        else:
            ax.set_ylabel("")
    
    # 隐藏多余的子图
    for i in range(n_interactions, rows * cols):
        row = i // cols
        col = i % cols
        axes[row, col].set_visible(False)
    
    # 添加总标题
    fig.suptitle("Band-Resolved pCOHP Analysis", fontsize=16, fontweight='bold', y=0.98)
    
    # 调整布局，为底部colorbar预留空间
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    
    # 添加水平colorbar在底部
    if scatter_plots:
        # 创建水平colorbar，放置在底部
        cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.03])  # [left, bottom, width, height]
        cbar = fig.colorbar(scatter_plots[0], cax=cbar_ax, orientation='horizontal')
        cbar.set_label("pCOHP (a.u.)", fontsize=14)
        cbar.ax.tick_params(labelsize=12)
    
    # 保存图片
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"绘图已保存为: {filename}")
    plt.show()

# --- 主程序 ---
if __name__ == "__main__":
    print("="*60)
    print("Band-Resolved pCOHP 分析脚本")
    print("="*60)
    
    # 1. 解析K点路径和高对称点
    print("\n步骤 1: 解析K点路径...")
    k_distances, high_sym_points = parse_kpoints_and_kpath()
    
    if k_distances is None:
        print("错误: 无法解析K点路径，程序终止。")
        exit()

    num_kpoints_from_kpoints = len(k_distances)
    print(f"✓ 成功解析 {num_kpoints_from_kpoints} 个K点")
    if high_sym_points:
        print(f"✓ 识别到 {len(high_sym_points)} 个高对称点")
    
    # 2. 查找并排序LOBSTER文件
    print("\n步骤 2: 查找LOBSTER文件...")
    file_list = sorted(glob.glob("KspaceCOHPBand*.lobster"), key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    if not file_list:
        print("错误: 在当前目录下没有找到任何 'KspaceCOHPBand*.lobster' 文件。")
        print("请确保LOBSTER输出文件在当前目录下。")
        exit()
    else:
        print(f"✓ 找到 {len(file_list)} 个能带文件")
        
        # 3. 处理LOBSTER文件并计算pCOHP
        print(f"\n步骤 3: 处理LOBSTER文件并计算pCOHP...")
        print(f"配置的相互作用类型: {len(interactions_to_plot)} 种")
        for interaction in interactions_to_plot:
            print(f"  - {interaction['label']}")
            
        all_pcohps_dict = {interaction['label']: [] for interaction in interactions_to_plot}
        all_energies = []
        labels = None
        
        for i, filename in enumerate(file_list):
            if (i + 1) % 10 == 0:
                print(f"  处理进度: {i+1}/{len(file_list)}")
                
            parsed_data = parse_single_band_cohp_file(filename)
            if parsed_data:
                num_kpoints_from_lobster, energies_one_band, cohps_one_band, current_labels = parsed_data
                
                if i == 0:
                    labels = current_labels
                    print(f"✓ 轨道标签数量: {len(labels)}")
                    if num_kpoints_from_lobster != num_kpoints_from_kpoints:
                        print("\n" + "="*50)
                        print("警告: K点数量不匹配!")
                        print(f"  - KPOINTS 文件中的路径点数: {num_kpoints_from_kpoints}")
                        print(f"  - LOBSTER 文件中的K点数: {num_kpoints_from_lobster}")
                        print("  - 将以 KPOINTS 文件中的点数为准进行截断。")
                        print("="*50 + "\n")

                if num_kpoints_from_lobster > num_kpoints_from_kpoints:
                    energies_one_band = energies_one_band[:num_kpoints_from_kpoints]
                    cohps_one_band = cohps_one_band[:num_kpoints_from_kpoints, :, :]
                elif num_kpoints_from_lobster < num_kpoints_from_kpoints:
                    print(f"错误: 文件 {filename} 的K点数少于KPOINTS文件，跳过。")
                    continue

                all_energies.append(energies_one_band)
                
                for interaction in interactions_to_plot:
                    pcohps_one_band = calculate_pcohp(cohps_one_band, labels, 
                                                      interaction["regex1"], 
                                                      interaction["regex2"])
                    all_pcohps_dict[interaction['label']].append(pcohps_one_band)
        
        if all_energies:
            print(f"✓ 成功处理 {len(all_energies)} 个能带")
            
            # 4. 保存数据和绘图
            print("\n步骤 4: 生成输出文件...")
            energies_arr = np.array(all_energies)
            pcohps_arr_dict = {label: np.array(pcohps) for label, pcohps in all_pcohps_dict.items()}
            
            combined_label = "_and_".join([inter['label'] for inter in interactions_to_plot])
            output_dat_filename = f"band_resolved_pcohp_{combined_label}.dat"
            output_png_filename = f"band_resolved_pcohp_{combined_label}.png"

            save_data_to_dat(output_dat_filename, k_distances, energies_arr, pcohps_arr_dict, fermi_energy=FERMI_ENERGY)
            
            plot_band_resolved_pcohp(output_png_filename, k_distances, energies_arr, pcohps_arr_dict, high_sym_points,
                                     fermi_energy=FERMI_ENERGY,
                                     energy_range=ENERGY_RANGE,
                                     vmax=VMAX)
            
            print(f"\n✓ 分析完成！输出文件:")
            print(f"  - 数据文件: {output_dat_filename}")
            print(f"  - 图片文件: {output_png_filename}")
            print("="*60)
        else:
            print("错误: 未能成功解析任何LOBSTER文件，程序终止。")
            print("请检查LOBSTER文件格式或轨道匹配设置。")

