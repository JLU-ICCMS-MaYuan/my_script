#!/usr/bin/env python3

import sys
import itertools
import collections
from typing import DefaultDict, Iterable, List, Tuple

import pandas as pd
import numpy as np
from pprint import pprint
from sympy import Matrix, lcm

from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry

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

def distinct_elems(mols: List[DefaultDict[str, int]]) -> List[str]:
    """Get the distinct elements that form the molecules."""
    elems = set()
    for mol in mols:
        for key in mol:
            elems.add(key)
    return list(elems)


def balance(left: List[DefaultDict[str, int]],
            right: List[DefaultDict[str, int]],
            verbose: bool) \
        -> Iterable[Tuple[np.ndarray, np.ndarray]]:
    """Balance the parsed left and ride sides."""
    elems = distinct_elems(left + right)
    if verbose:
        print('Distinct elements:', ', '.join(elems))
    lin_sys = np.zeros((len(elems), len(left) + len(right)), dtype=int)
    for idx_elem, elem in enumerate(elems):
        for idx_mol, mol in enumerate(itertools.chain(left, right)):
            lin_sys[idx_elem, idx_mol] = mol[elem]
    lin_sys[:, len(left):] *= -1
    if verbose:
        print('Linear system of equations matrix', lin_sys, sep='\n')
    kernel = Matrix(lin_sys).nullspace(simplify=True)
    if verbose:
        print('Nullity =', len(kernel))
    for ker in kernel:
        numers = np.array([num.as_numer_denom()[0] for num in ker])
        denoms = np.array([num.as_numer_denom()[1] for num in ker])
        lcm = np.lcm.reduce(denoms)
        coefs = numers * lcm / denoms
        coefs /= np.gcd.reduce(coefs)
        if np.any(coefs < 0) and np.any(coefs > 0):
            continue
        coefs = np.abs(coefs)
        if verbose:
            print('Kernel basis', ker, 'simplified to', coefs)
    yield coefs[:len(left)], coefs[len(left):]

def get_chemical_equation(coefs_left, entries):
    """Construct string of coefficients and molecules."""
    multipled = []
    # print(coefs_left, formula); input()
    for coef, entry in zip(coefs_left, entries):
        formula = entry.composition.formula.replace(' ', '')
        entryid = entry.entry_id
        speciel_formula = 'id' + str(entryid) + '_' + formula
        multipled.append(f'{coef}*{speciel_formula}' if coef != 1 else speciel_formula)
    return ' + '.join(multipled)

def compute_form_energy(left_coef, left_compound, right_coef, right_compound):
    # product_idx 是生成物在product_reacants_entries中的索引号组成的列表
    form_energy = 0
    left_total_numatm = 0
    right_total_numatm = 0
    # print(left_compound, right_compound); input()
    for coef, mol in zip(left_coef, left_compound):
        # print(coef, mol); input()
        natm = mol.composition.num_atoms
        energy = mol.energy
        left_total_numatm += natm*coef
        form_energy  += natm*coef*energy

    for coef, mol in zip(right_coef, right_compound):
        natm = mol.composition.num_atoms
        energy = mol.energy
        right_total_numatm += natm*coef
        form_energy  += -natm*coef*energy
    
    if left_total_numatm == right_total_numatm and left_total_numatm != 0:
        form_energy_peratom = 1000*form_energy/right_total_numatm
    elif left_total_numatm == right_total_numatm and left_total_numatm == 0:
        form_energy_peratom = None
    else:
        print(f"The left part numberofatom:{left_total_numatm} != the left part numberofatom:{right_total_numatm}")
        sys.exit
    return form_energy_peratom



if __name__ == "__main__":
    # 如果生成物有n个，就放在列表的最后n列
    enthalpy_csvfile =     sys.argv[1]
    numberofproducts = int(sys.argv[2])
    numberofreactions= int(sys.argv[3])
    enthalpy_datas   = pd.read_csv(enthalpy_csvfile, index_col=0, header=1, na_values=['--'])
    result_dt = collections.defaultdict(list)
    presses = []
    for press, idxcomp_energy in enthalpy_datas.iterrows():
        presses.append(press)
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
                formula1 = compound1.composition.formula.replace(' ', '')
                formula2 = compound2.composition.formula.replace(' ', '')
                if formula1 == formula2:
                    # print(f'formula1 = {formula1}, its entry_id = {compound1.entry_id}')
                    # print(f'formula2 = {formula2}, its entry_id = {compound2.entry_id}')
                    break
            else:
                combination_reacants.append(list(combination))
        
        # 拿到所有反应物的排列组合
        for reactions in combination_reacants:
            left_compound  = [comp.composition.get_el_amt_dict() for comp in reactions]
            right_compound = [comp.composition.get_el_amt_dict() for comp in product]
            left_formula   = [comp.composition.formula.replace(' ', '') for comp in reactions]
            right_formula  = [comp.composition.formula.replace(' ', '') for comp in product]
            solutions = balance(left_compound, right_compound, verbose=False)
            for left_coef, right_coef in solutions:
                left_part  = get_chemical_equation(left_coef, reactions)
                right_part = get_chemical_equation(right_coef, product)
                chemical_equation = left_part + ' = ' + right_part
                # print(left_part, '->', right_part); input()
                form_energy_peratom = compute_form_energy(left_coef, reactions, right_coef, product)
                if form_energy_peratom is not None:
                    result_dt[chemical_equation].append(form_energy_peratom)
    
    # print(result_dt)
    total_result_pd = pd.DataFrame(data=result_dt, index=presses)
    total_result_pd.to_csv("formed-enthalpy.csv")
    




