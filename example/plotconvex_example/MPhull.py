from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.analysis.phase_diagram import PDPlotter
from pymatgen.ext.matproj import MPRester
from IPython.core.display import display
import pandas as pd
import re


def get_above_hull_and_formation_enrgy_by_total_energy(tot_energy_excel_path, MPAPI_KEY):
    import pandas as pd
    # 导入化学式——总能量数据表
    tot_energy_excel = pd.read_excel(tot_energy_excel_path)
    formu_energy_dic = dict(zip(tot_energy_excel.iloc[:, 0], tot_energy_excel.iloc[:, 1]))
    # 批量处理分支
    # 获取MP数据
    pd_list = []
    forma_energy_list = []
    e_above_hull_list = []
    entries = []
    with MPRester(MPAPI_KEY) as mpr:
        for key in formu_energy_dic:
            div = re.compile(r"([A-Z][a-z])|([A-Z])|(\(.*\))")
            ele_list = div.findall(key)
            formu_key = ''
            for el in ele_list:
                for e in el:
                    if e != '':
                        formu_key += e + '-'
            formu_key = formu_key.strip('-')
            print(key)
            #entries = mpr.get_entries_in_chemsys(formu_key)  # 获取MP数据库所有相关相
            pd_entry = PDEntry(composition=key, energy=formu_energy_dic[key])  # 创建自己导入的化学式-总能相
            entries.append(pd_entry)
            # 画相图
            phase_diagram = PhaseDiagram(entries=entries)

            # 计算导入相形成能
            forma_energy_list.append(phase_diagram.get_form_energy_per_atom(pd_entry))
            # 计算导入相above hull
            e_above_hull_list.append(phase_diagram.get_e_above_hull(pd_entry))

            pd_list.append(phase_diagram)
    tot_energy_excel['formation_energy_per_atom'] = forma_energy_list
    tot_energy_excel['e_above_hull'] = e_above_hull_list
    res_excel = tot_energy_excel
    display(res_excel)
    return pd_list, formu_energy_dic, res_excel

if __name__ == '__main__':
    pd_list, formu_energy_dic, res_excel = get_above_hull_and_formation_enrgy_by_total_energy(tot_energy_excel_path='ehull.xlsx', MPAPI_KEY='LmfQdNl7YVNPXBG1ncy8OJqlkOMylntj')

    # 将结果保存为Excel
    res_excel.to_excel('aboveHull_result.xlsx')
    
    # 画相图，只能手动更改index值指定画第几个化学式的相图
    #index = 1  # 第index个样本的相图
    #display(res_excel.iloc[index, :])
    #plotter = PDPlotter(pd_list[index])
    #plotter.get_plot()
