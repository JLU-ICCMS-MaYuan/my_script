#!/usr/bin/env python3
"""
LOBSTER轨道相互作用分析工具
读取ICOHPLIST.lobster, ICOOPLIST.lobster, ICOBILIST.lobster文件
根据用户指定的轨道对计算相互作用值
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import defaultdict
from pathlib import Path


class LobsterAnalyzer:
    """LOBSTER数据分析器"""
    
    def __init__(self, lobster_dir="."):
        """
        初始化分析器
        
        Args:
            lobster_dir: LOBSTER文件所在目录
        """
        self.lobster_dir = Path(lobster_dir)
        self.data = {}
        
    def read_lobster_files(self):
        """读取三个LOBSTER文件"""
        file_types = {
            'icohp': 'ICOHPLIST.lobster',
            'icoop': 'ICOOPLIST.lobster', 
            'icobi': 'ICOBILIST.lobster'
        }
        
        for file_type, filename in file_types.items():
            filepath = self.lobster_dir / filename
            if not filepath.exists():
                raise FileNotFoundError(f"找不到文件: {filepath}")
            
            print(f"正在读取 {filename}...")
            self.data[file_type] = self._parse_lobster_file(filepath)
            
    def _parse_lobster_file(self, filepath):
        """解析单个LOBSTER文件"""
        data = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        # 跳过前两行标题
        for i, line in enumerate(lines[2:], start=3):
            if line.strip() == "":
                continue
                
            parts = line.strip().split()
            if len(parts) < 6:
                continue
                
            # 解析数据行
            cohp_num = int(parts[0])
            atom_mu = parts[1]
            atom_nu = parts[2] 
            distance = float(parts[3])
            translation = tuple(map(int, parts[4:7]))
            value = float(parts[7])
            
            data.append({
                'cohp_num': cohp_num,
                'atom_mu': atom_mu,
                'atom_nu': atom_nu,
                'distance': distance,
                'translation': translation,
                'value': value
            })
            
        return pd.DataFrame(data)
    
    def expand_orbital_spec(self, orbital_list, data_type='icohp'):
        """
        展开轨道规范，支持通配符匹配
        
        Args:
            orbital_list: 轨道列表，如["La_4f", "H1_1s"]
            data_type: 数据类型
            
        Returns:
            匹配的完整轨道列表
        """
        df = self.data[data_type]
        all_orbitals = set(df['atom_mu'].tolist() + df['atom_nu'].tolist())
        
        expanded = []
        for orbital_spec in orbital_list:
            matched = self._match_orbital_pattern(orbital_spec, all_orbitals)
            expanded.extend(matched)
            
        return list(set(expanded))  # 去重
    
    def _match_orbital_pattern(self, pattern, all_orbitals):
        """匹配轨道模式"""
        matched = []
        
        if pattern in all_orbitals:
            # 精确匹配
            return [pattern]
            
        # 通配符匹配
        if '_' in pattern:
            # 如 "La_4f" 匹配所有La原子的4f轨道
            parts = pattern.split('_')
            element = parts[0]
            orbital_type = '_'.join(parts[1:]) if len(parts) > 1 else ""
            
            for orbital in all_orbitals:
                if orbital.startswith(element):
                    if orbital_type == "":
                        matched.append(orbital)
                    elif orbital_type in orbital:
                        matched.append(orbital)
        else:
            # 如 "H" 匹配所有H原子的轨道
            for orbital in all_orbitals:
                if orbital.startswith(pattern):
                    matched.append(orbital)
                    
        return matched
    
    def calculate_interaction(self, orbital_pairs, data_type='icohp'):
        """
        计算轨道对相互作用
        
        Args:
            orbital_pairs: 轨道对规范，如[["La1_5p_y", "La1_5p_z"], ["H1_1s", "H2_1s"]]
            data_type: 数据类型
            
        Returns:
            计算结果和标签
        """
        left_orbitals, right_orbitals = orbital_pairs
        
        # 展开轨道规范
        left_expanded = self.expand_orbital_spec(left_orbitals, data_type)
        right_expanded = self.expand_orbital_spec(right_orbitals, data_type)
        
        print(f"左侧轨道: {left_expanded}")
        print(f"右侧轨道: {right_expanded}")
        
        df = self.data[data_type]
        
        # 收集所有匹配的轨道对及其键编号
        matched_bonds = set()  # 用于统计涉及的键数量
        total_value = 0
        
        for left_orbital in left_expanded:
            for right_orbital in right_expanded:
                # 查找相互作用值（双向查找）
                mask1 = (df['atom_mu'] == left_orbital) &  (df['atom_nu'] == right_orbital)
                mask2 = (df['atom_mu'] == right_orbital) & (df['atom_nu'] == left_orbital)
                
                matched_rows = df[mask1 | mask2]
                if not matched_rows.empty:
                    # 累加所有匹配的值
                    total_value += matched_rows['value'].sum()
                    # 记录涉及的键编号
                    matched_bonds.update(matched_rows['cohp_num'].tolist())
        
        print(f"匹配的轨道对总和: {total_value}")
        print(f"涉及的键数量: {len(matched_bonds)}")
        print(f"涉及的键编号: {sorted(matched_bonds)}")
        
        # 按键的数量求平均
        if len(matched_bonds) > 0:
            result = total_value / len(matched_bonds)
        else:
            result = 0.0
            
        # 生成标签
        label = self._generate_label(orbital_pairs)
        
        return result, label
    
    def _group_by_atoms(self, orbitals):
        """按原子分组轨道"""
        atoms = defaultdict(list)
        for orbital in orbitals:
            # 提取原子名称（如La1, H2等）
            atom_match = re.match(r'([A-Za-z]+\d*)', orbital)
            if atom_match:
                atom = atom_match.group(1)
                atoms[atom].append(orbital)
        return dict(atoms)
    
    def _generate_label(self, orbital_pairs):
        """生成输出标签"""
        left_orbitals, right_orbitals = orbital_pairs
        
        # 简化轨道描述
        left_label = self._simplify_orbital_list(left_orbitals)
        right_label = self._simplify_orbital_list(right_orbitals)
        
        return f"{left_label}_{right_label}"
    
    def _simplify_orbital_list(self, orbital_list):
        """简化轨道列表为标签"""
        if len(orbital_list) == 1:
            orbital = orbital_list[0]
            if '_' in orbital:
                return orbital.replace('_', '_').replace('^', '^')
            return orbital
        
        # 多个轨道的情况，尝试找共同模式
        common_prefix = self._find_common_prefix(orbital_list)
        if common_prefix:
            return common_prefix
        
        # 如果没有共同前缀，连接所有轨道
        return '_'.join(orbital_list)
    
    def _find_common_prefix(self, orbital_list):
        """找到轨道列表的共同前缀"""
        if not orbital_list:
            return ""
            
        # 尝试找元素和轨道类型的共同模式
        elements = set()
        orbital_types = set()
        
        for orbital in orbital_list:
            parts = orbital.split('_')
            if len(parts) >= 2:
                # 提取原子（如La1）和轨道类型（如4f）
                atom_match = re.match(r'([A-Za-z]+)', parts[0])
                if atom_match:
                    elements.add(atom_match.group(1))
                    
                orbital_type_match = re.match(r'(\d+[spdf])', parts[1])
                if orbital_type_match:
                    orbital_types.add(orbital_type_match.group(1))
        
        if len(elements) == 1 and len(orbital_types) == 1:
            return f"{list(elements)[0]}_{list(orbital_types)[0]}"
        elif len(elements) == 1:
            return list(elements)[0]
            
        return ""
    
    def process_multiple_pairs(self, orbital_pairs_list):
        """处理多个轨道对"""
        results = {
            'icohp': {},
            'icoop': {},
            'icobi': {}
        }
        
        for orbital_pairs in orbital_pairs_list:
            print(f"\n处理轨道对: {orbital_pairs}")
            
            for data_type in ['icohp', 'icoop', 'icobi']:
                value, label = self.calculate_interaction(orbital_pairs, data_type)
                results[data_type][label] = value
                print(f"{data_type.upper()}: {label} = {value:.6f}")
                
        return results
    
    def save_results(self, results):
        """保存结果到文件"""
        for data_type, data in results.items():
            filename = f"proj_{data_type}.dat"
            filepath = self.lobster_dir / filename
            
            with open(filepath, 'w') as f:
                for label, value in data.items():
                    f.write(f"{label}\t{value:.6f}\n")
            
            print(f"已保存 {filename}")
    
    def plot_results(self, results):
        """绘制结果图"""
        labels = list(results['icohp'].keys())
        icohp_values = [results['icohp'][label] for label in labels]
        icoop_values = [results['icoop'][label] for label in labels]
        icobi_values = [results['icobi'][label] for label in labels]
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        # ICOHP
        bars1 = ax1.bar(range(len(labels)), icohp_values)
        ax1.set_xlabel('Orbital Pairs')
        ax1.set_ylabel('ICOHP')
        ax1.set_xticks(range(len(labels)))
        ax1.set_xticklabels(labels, rotation=45, ha='right')
        
        # 在柱子上显示数值
        for bar, value in zip(bars1, icohp_values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.4f}', ha='center', va='bottom' if height > 0 else 'top')
        
        # ICOOP
        bars2 = ax2.bar(range(len(labels)), icoop_values)
        ax2.set_xlabel('Orbital Pairs')
        ax2.set_ylabel('ICOOP')
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=45, ha='right')
        
        # 在柱子上显示数值
        for bar, value in zip(bars2, icoop_values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.5f}', ha='center', va='bottom' if height >= 0 else 'top')
        
        # ICOBI
        bars3 = ax3.bar(range(len(labels)), icobi_values)
        ax3.set_xlabel('Orbital Pairs')
        ax3.set_ylabel('ICOBI')
        ax3.set_xticks(range(len(labels)))
        ax3.set_xticklabels(labels, rotation=45, ha='right')
        
        # 在柱子上显示数值
        for bar, value in zip(bars3, icobi_values):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.4f}', ha='center', va='bottom' if height >= 0 else 'top')
        
        plt.tight_layout()
        
        plot_path = self.lobster_dir / "orbital_interactions.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"已保存图像: {plot_path}")


def main():
    """主函数"""
    # 示例用法
    analyzer = LobsterAnalyzer(".")
    
    # 读取LOBSTER文件
    analyzer.read_lobster_files()
    
    # 定义轨道对 - 测试具体指定轨道
    
    # orbital_pairs_list = [
    #         [["La_4f"], ["H_1s"]],  # 测试有值的轨道对
    #     ]
    
    # orbital_pairs_list = [
    #         [["La1_4f_y(3x^2-y^2)"], ["H_1s"]],
    #         [["La1_4f_xyz"],         ["H_1s"]],
    #         [["La1_4f_yz^2"],        ["H_1s"]],
    #         [["La1_4f_z^3"],         ["H_1s"]],
    #         [["La1_4f_xz^2"],        ["H_1s"]],
    #         [["La1_4f_z(x^2-y^2)"],  ["H_1s"]],
    #         [["La1_4f_x(x^2-3y^2)"], ["H_1s"]],
    #     ]
    
    # orbital_pairs_list = [
    #     [["La1_4f_y(3x^2-y^2)"], ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_xyz"],         ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_yz^2"],        ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_z^3"],         ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_xz^2"],        ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_z(x^2-y^2)"],  ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    #     [["La1_4f_x(x^2-3y^2)"], ["H4_1s", "H5_1s", "H6_1s", "H7_1s", "H8_1s", "H9_1s", "H16_1s", "H17_1s", "H18_1s", "H19_1s", "H20_1s", "H21_1s"]],
    # ]
    
    orbital_pairs_list = [
        [["La1_4f_y(3x^2-y^2)"], ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_xyz"],         ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_yz^2"],        ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_z^3"],         ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_xz^2"],        ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_z(x^2-y^2)"],  ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
        [["La1_4f_x(x^2-3y^2)"], ["H10_1s", "H11_1s", "H12_1s", "H13_1s", "H14_1s", "H15_1s"]],
    ]
    
    # orbital_pairs_list = [
    #     [["La1_4f_y(3x^2-y^2)"], ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_xyz"],         ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_yz^2"],        ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_z^3"],         ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_xz^2"],        ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_z(x^2-y^2)"],  ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    #     [["La1_4f_x(x^2-3y^2)"], ["H22_1s", "H23_1s", "H24_1s", "H25_1s", "H26_1s", "H27_1s"]],
    # ]
    
    # 处理轨道对
    results = analyzer.process_multiple_pairs(orbital_pairs_list)
    
    # 保存结果
    analyzer.save_results(results)
    
    # 绘制图像
    analyzer.plot_results(results)


if __name__ == "__main__":
    main()