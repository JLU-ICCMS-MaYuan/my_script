#!/usr/bin/env python3
import os
import numpy as np

formulas = np.genfromtxt('formula_pool', dtype=int)
print('Before filter, number of formulas: {}'.format(formulas.shape[0]))


# 备份读入的formula_pool，因为最后生成文件要覆盖formula_pool的内容
i=1; formulas_backup = 'formula_pool_backup'+str(i)
while os.path.exists(formulas_backup):
    i+=1
    formulas_backup = 'formula_pool_backup'+str(i)
np.savetxt(formulas_backup, formulas, fmt='%d')


# 筛选目标formulas
filtered_formulas = formulas[(formulas[:, 0] <= 6) & (formulas[:, 1] <= 6)]
print('After filter, number of formulas: {}'.format(filtered_formulas.shape[0]))
# 扩胞目标formulas
extended_formulas = np.vstack([filtered_formulas, filtered_formulas*2, filtered_formulas*3, filtered_formulas*4])
print('After vary cell, number of formulas: {}'.format(extended_formulas.shape[0]))
# 整合所有formula
total_natoms = np.sum(extended_formulas, axis=1)
print('max total number of atoms: {}'.format(np.max(total_natoms)))


# 覆盖formula_pool
np.savetxt('formula_pool', extended_formulas, fmt='%d')