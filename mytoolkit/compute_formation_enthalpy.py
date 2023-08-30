#!/usr/bin/env python3

import sys
import itertools
import collections
import time

import pandas as pd
import numpy as np
from pprint import pprint
from sympy import Matrix, lcm, Rational

from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram

# 如果你要使用这个脚本，你需要按照如下的格式准备数据，
# 1. 这个数据格式必须是csv格式，
# 2. 第一行是你自己定义的各个结构的名称，随便怎么写都可以，因为这个脚本不会读入第一行
# 3. 第二行是脚本要求你重新定义的各个结构的名称，必须为`数字-化学式`，数字必须是列数，例如1.H-P63mc在第一列就是1-H
# 4. 你可以不按照如下对的这么整齐，因为csv默认按照逗号分隔，不是空格分割。
# 5. 最后几列必须是生成物
# 6. 输出的文件叫formed-enthalpy.csv
# 7. 输入参数有3个，第1个是你准备的焓值文件，第2个是指定生成物的数量，一个是指定反应物的数量
"""
name,       1.H-P63mc,   21.LaH3,    24.YH2-Fm-3m, 18.CeH3-50GPa, 33.ThH4-15-85-P321, 001-BeH2,      LaYCeThBe4H32
Press(GPa), 1-H,         3-LaH3,     5-YH2,        8-CeH3,        14-ThH4,            19-BeH2,       20-LaYCeThBe4H32
10,        -2.95760473, -3.4877319, -5.13175351,  -3.92713112,   -4.05826962,        -5.51291338,   -3.95950239                                                                                                       
20,        -2.6713429,  -2.93698821,-4.49373914,  -3.40366919,   -3.57765997,        -5.20412224,   -3.63861714
30,        -2.4311295,  -2.42867868,-3.89668335,  -2.91361873,   -3.12443305,        -4.92204357,   -3.3365023
40,        -2.21776141, -1.95137263,-3.33129567,  -2.44992903,   -2.69310992,        -4.65846137,   -3.04918475
50,        -2.02286976, -1.49843756,-2.79161411,  -2.00726826,   -2.28001663,        -4.40929342,   -2.77403504
60,        -1.84182373, -1.06551571,-2.27349076,  -1.58248683,   -1.88251391,        -4.17179548,   -2.50914441
70,        -1.67171379, -0.64952649,-1.77388899,  -1.17315931,   -1.49860812,        -3.94408329,   -2.25310365
80,        -1.51055794, -0.24817143,-1.29049139,  -0.77740955,   -1.12674552,        -3.72478002,   -2.00483278
90,        -1.35692834,  0.14029836,-0.82150196,  -0.39375822,   -0.76567542,        -3.51283738,   -1.76347291
100,       -1.20975282,  0.51728745,-0.36542815,  -0.0210147,    -0.41436475,        -3.30741971,   -1.52832805
110,       -1.06822025,  0.88423433, 0.07933491,   0.3417946,    -0.07195121,         --,           -1.29880992
120,       -0.93164902,  1.24124169, 0.51308954,   0.69550899,    0.2623178,         -2.91307268,   -1.0744323
130,       -0.79950153,  1.58961144, 0.93698481,   1.04083835,    0.58907071,         --,           -0.85480325
140,       -0.67133393,  1.93011396, 1.35174782,   1.37840671,    0.90886508,         --,           -0.63954283
150,       -0.54677398,  2.26341065, 1.75800913,   1.70864529,    1.22218614,         -2.35773708,  -0.42832854
160,       -0.42550605,  2.59001417, 2.15630701,   2.03158367,    1.52946701,         --,           -0.22089103
170,       -0.30725955,  2.9113808,  2.54663092,   2.34902635,    1.83108967,         --,           -0.01698845
180,       -0.1917999,   3.2249112,  2.92936596,   2.65994143,    2.1273921,          -1.83661578,   0.18359596
190,       -0.07892249,  3.53398105, 3.30520196,   2.96494061,    2.418685,           --,            0.38144126
200,        0.0315527,   3.83767844, 3.67441807,   3.26454873,    2.70524702,         -1.50507977,   0.57597564
"""

def make_elementMatrix(product_reacants_entries):
    # input(product_reacants_entries)
    elementmatrix = []
    elementlist   = []
    for idx, entry in enumerate(product_reacants_entries):
        if(idx == len(elementmatrix)):
            elementmatrix.append([])
        for x in elementlist:
            elementmatrix[idx].append(0)
        # print(f"第一次补零: elementmatrix={elementmatrix}"); input()
        el_amt_dict = entry.composition.get_el_amt_dict()
        for ele, amt in el_amt_dict.items():
            if (ele not in elementlist):
                elementlist.append(ele)
                for i in range(len(elementmatrix)):
                    elementmatrix[i].append(0)
            column=elementlist.index(ele)
            elementmatrix[idx][column]+=int(amt)
        # pprint(entry)
        # print(f"第2次: elementmatrix={elementmatrix}"); input()
    return elementmatrix

def get_chemical_equation(coeffients, product_reacants_entries, number_of_products):
    product_entries   = product_reacants_entries[-number_of_products:]
    product_coeffient = coeffients[-number_of_products:]
    reactant_entries  = product_reacants_entries[:-number_of_products]
    reactant_coeffient= coeffients[:-number_of_products]

    reactants_items = []
    for rea_coef, rea_entry in zip(reactant_coeffient, reactant_entries):
        formula = rea_entry.composition.formula.replace(' ', '')
        entryid = rea_entry.entry_id
        reactants_items.append(str(int(rea_coef)) +'*'+ '('+str(entryid)+'-'+formula+')')
    reactants_strings ='+ '.join(reactants_items)

    products_items = []
    for rea_coef, prod_entry in zip(product_coeffient, product_entries):
        formula = prod_entry.composition.formula.replace(' ', '')
        entryid = prod_entry.entry_id
        products_items.append(str(int(rea_coef)) +'*'+ '('+str(entryid)+'-'+formula+')')
    products_strings ='+ '.join(products_items)

    chemical_equation = reactants_strings + " -> " + products_strings
    # print(chemical_equation)
    return chemical_equation

def compute_form_energy(coeffients, product_reacants_entries, number_of_products):
    # product_idx 是生成物在product_reacants_entries中的索引号组成的列表
    product_idx = list(range((len(product_reacants_entries)-number_of_products), len(product_reacants_entries)))
    # print(product_idx); input("product_idx") # 输出生成物的索引号
    # product_numatm = 0
    # reactant_numatm= 0
    # for prod_entry in product_reacants_entries[-number_of_products:]:
    #     natm = prod_entry.composition.num_atoms
    #     product_numatm += natm
    # for react_entry in product_reacants_entries[:-number_of_products]:
    #     natm = react_entry.composition.num_atoms
    #     reactant_numatm += natm
    # if product_numatm == reactant_numatm:
    #     total_numatm = product_numatm
    # else:
    #     print(f"product_numatm={product_numatm}")
    #     print(f"reactant_numatm={reactant_numatm}")
    #     print("Number of atoms in product and reactant are not equal. I will exit !!!")
    #     # sys.exit(1)
    form_energy = 0
    total_numatm = 0
    for prod_entry in product_reacants_entries[-number_of_products:]:
        natm = prod_entry.composition.num_atoms
        total_numatm += natm
    for idx, prod_reac in enumerate(product_reacants_entries):
        numatm = prod_reac.composition.num_atoms
        energy = prod_reac.energy
        if idx not in product_idx :
            # print(f"numatm={numatm}, coeffients[{idx}]={coeffients[idx]}")
            form_energy += numatm*coeffients[idx]*energy
        else:
            form_energy += -numatm*coeffients[idx]*energy
    form_energy_peratom = 1000*form_energy/total_numatm
    # print(f"form_energy_peratom={form_energy_peratom}")#; input()
    return form_energy_peratom



if __name__ == "__main__":
    # 如果生成物有n个，就放在列表的最后n列
    enthalpy_csvfile =     sys.argv[1]
    numberofproducts = int(sys.argv[2])
    numberofreactions= int(sys.argv[3])
    enthalpy_datas   = pd.read_csv(enthalpy_csvfile, index_col=0, header=1, na_values=['--'])
    result_dt = collections.defaultdict(list)
    for press, idxcomp_energy in enthalpy_datas.iterrows():
        print(f"press={press}")
        ini_entries = []
        for idxcomp, enthalpy in idxcomp_energy.items():
            if enthalpy != np.nan:
                _comp  = Composition(idxcomp.split("-")[-1])
                _entry = PDEntry(_comp, enthalpy)
                _entry.entry_id = idxcomp.split("-")[0]
                ini_entries.append(_entry)
            else:
                print(f"{press}GPa下部分结构未计算能量")
                break

        product    = ini_entries[-numberofproducts:]
        reactants  = ini_entries[:-numberofproducts]
        combination_reacants = []
        for combination in list(itertools.combinations(reactants, numberofreactions)):
            for compound1, compound2 in list(itertools.combinations(combination, 2)):
                if compound1.composition.elements == compound2.composition.elements:
                    break
            else:
                combination_reacants.append(list(combination))
        
        # 拿到所有反应物的排列组合
        prods_reacs = [rec + product for rec in combination_reacants] 
        for prod_reac in prods_reacs:
            # print(prod_reac); input()
            ele_matrix = make_elementMatrix(prod_reac)
            ele_matrix = Matrix(ele_matrix)
            ele_matrix = ele_matrix.transpose()
            solution   = ele_matrix.nullspace()[0]
            multiple   = lcm([val.q for val in solution]) # 获得所有系数的分母的值，然后取它们的最小公倍数
            coeffients = abs(multiple*solution)
            # print(coeffients); input()
            chemical_equation   = get_chemical_equation(coeffients, prod_reac, numberofproducts)
            # print(chemical_equation); input()
            form_energy_peratom = compute_form_energy(coeffients, prod_reac, numberofproducts)
            # print(form_energy_peratom); input()
            # parital_result_pd.at[press, chemical_equation] = form_energy_peratom
            result_dt[chemical_equation].append(form_energy_peratom)
            # print(parital_result_pd);input()
        # total_result_pd = pd.concat([total_result_pd, parital_result_pd], axis=1)
    
    total_result_pd = pd.DataFrame(data=result_dt, index=[p for p in range(10, 210, 10)]+[250, 300])
    total_result_pd.to_csv("formed-enthalpy.csv")
    




