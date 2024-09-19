#!/usr/bin/env python3
import pandas as pd
import re
import sys


def read_compounds_from_file(file_name):
    # 使用 pandas 读取文件
    df = pd.read_table(file_name, header=None, names=['formula', 'min', 'max'], sep='\s+')

    compounds = {}
    ranges = {}
    
    for _, row in df.iterrows():
        compound = row['formula']  # 例如 CeH10
        start_range = int(row['min'])   # 最小数量
        end_range = int(row['max'])     # 最大数量
        # 使用正则表达式拆解化合物，例如 CeH10 -> {'Ce': 1, 'H': 10}
        elements = {}
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', compound)
        # [A-Z] 匹配元素的第一个大写字母
        # [a-z]* 匹配元素名称中可能存在的小写字母
        # (\d*) 匹配可能存在的数字（即原子数），如果没有数字，则默认为1
        for element, count in matches:
            elements[element] = int(count) if count else 1

        # 存储化合物及其数量范围
        compounds[compound] = elements
        ranges[compound] = range(start_range, end_range+1)
    
    return compounds, ranges


if __name__ == "__main__":
    info = '''You can't write `binary_pools.dat` like:
First Type:
Ce1H10 1 3
Sc1H3 1 3
H1  0  8

Second Type:
Ce1H10 0 3
Ce1H9 0 3
Sc1H3 0 3
Sc1H6 0 3

Second Type:
Ce1H10 0 3
Ce1H9 0 3
Sc1H3 0 3
Sc1H6 0 3
H  0 3

Then you can run it by: getTernaryFromCombinateBinary.py Ce Sc H
Finally, it will create two file: ` formula_pool_for_magus` and `compound_combinations.csv`
'''
    print(info)
    try:
        prefered_order = sys.argv[1:4]
    except:
        prefered_order = []
        print("You didn't specify prefered order of elements. It will use default orders")
    # 从文件读取化合物数据
    compounds, ranges = read_compounds_from_file('binary_pools.dat')
    # print(compounds) # {'Ce1H10': {'Ce': 1, 'H': 10}, 'Sc1H3': {'Sc': 1, 'H': 3}, 'H1': {'H': 1}}
    # print(ranges)    # {'Ce1H10': range(1, 4), 'Sc1H3': range(1, 4), 'H1': range(0, 9)}
    # 生成配比组合
    from itertools import product

    # 将化合物和范围进行组合
    compound_names = list(compounds.keys()) # 获得现有的所有化学式：CeH10, ScH3, H
    compound_combinations = product(*(ranges[name] for name in compound_names)) # 获得每个化学式对应的组分变化范围: range(1, 4), range(1, 4), range(0, 9)
    # 关键点来了：
    #           (ranges[name] for name in compound_names)  =         (range(1, 4), range(1, 4), range(0, 9))
    #          *(ranges[name] for name in compound_names)  =          range(1, 4), range(1, 4), range(0, 9)   # 星号的作用是解包，就是去除括号
    #  product(*(ranges[name] for name in compound_names)) = product( range(1, 4), range(1, 4), range(0, 9) )
    
    # 打开文件并写入数据
    numbers = []
    with open('compound_combinations.csv', 'w') as f:
        # 写入第一行，指明二元化合物的名称
        f.write(','.join(compound_names) + ',Formula\n')

        for combination in compound_combinations:
            # print(combination) # (2, 2, 6)
            # 初始化新化合物的组成
            new_compound = {}

            # 遍历每个化合物及其对应的数量
            for i, compound_name in enumerate(compound_names):
                compound = compounds[compound_name] # compound = CeH10
                count = combination[i]  # count = 2

                # 计算新化合物的元素数量
                for element, element_count in compound.items(): # element=Ce element_count=10
                    if element in new_compound:
                        new_compound[element] += count * element_count
                    else:
                        new_compound[element] = count * element_count

            if 0 in list(new_compound.values()):
                continue

            # 生成化学式
            if prefered_order:
                formula = ''.join([f'{element}{new_compound[element]}' for element in prefered_order if element in new_compound])
                number  = '    '.join([f'{new_compound[element]}' for element in prefered_order if element in new_compound])
            else:
                formula = ''.join([f'{element}{new_compound[element]}' for element in new_compound])
                number = '    '.join([f'{new_compound[element]}' for element in new_compound])
            
            numbers.append(number)
            
            # 输出结果
            # print(f'{combination} -> {formula}')
            # 写入文件：组合值 + 化学式
            f.write(','.join(map(str, combination)) + f',{formula}\n')

    with open('formula_pool_for_magus', 'w') as f:
        for number in numbers:
            f.write(f'{number}\n')