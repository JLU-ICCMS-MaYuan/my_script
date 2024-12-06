# 字典：元素对应的赝势文件路径
base_pseudo_path = "/home/jildxwlxyljlstdui/cccs-share02/lijx/liz/potpaw_PBE/"
pseudo_paths = {
    'H': 'H/POTCAR',
    'Li': 'Li_sv/POTCAR',
    'Be': 'Be_sv/POTCAR',
    'B': 'B/POTCAR',
    'C': 'C/POTCAR',
    'N': 'N/POTCAR',
    'O': 'O/POTCAR',
    'F': 'F/POTCAR',
    'Ne': 'Ne/POTCAR',
    'Na': 'Na_pv/POTCAR',
    'Mg': 'Mg_pv/POTCAR',
    'Al': 'Al/POTCAR',
    'Si': 'Si/POTCAR',
    'P': 'P/POTCAR',
    'S': 'S/POTCAR',
    'Cl': 'Cl/POTCAR',
    'K': 'K_sv/POTCAR',
    'Ca': 'Ca_sv/POTCAR',
    'Sc': 'Sc_sv/POTCAR',
    'Ti': 'Ti_sv/POTCAR',
    'V': 'V_sv/POTCAR',
    'Cr': 'Cr_pv/POTCAR',
    'Mn': 'Mn_pv/POTCAR',
    'Fe': 'Fe_pv/POTCAR',
    'Co': 'Co/POTCAR',
    'Ni': 'Ni_pv/POTCAR',
    'Cu': 'Cu_pv/POTCAR',
    'Zn': 'Zn/POTCAR',
    'Ga': 'Ga_d/POTCAR',
    'Ge': 'Ge_d/POTCAR',
    'As': 'As/POTCAR',
    'Se': 'Se/POTCAR',
    'Br': 'Br/POTCAR',
    'Kr': 'Kr/POTCAR',
    'Rb': 'Rb_sv/POTCAR',
    'Sr': 'Sr_sv/POTCAR',
    'Y': 'Y_sv/POTCAR',
    'Zr': 'Zr_sv/POTCAR',
    'Nb': 'Nb_sv/POTCAR',
    'Mo': 'Mo_pv/POTCAR',
    'Tc': 'Tc_pv/POTCAR',
    'Ru': 'Ru_pv/POTCAR',
    'Rh': 'Rh_pv/POTCAR',
    'Pd': 'Pd/POTCAR',
    'Ag': 'Ag/POTCAR',
    'Cd': 'Cd/POTCAR',
    'In': 'In_d/POTCAR',
    'Sn': 'Sn_d/POTCAR',
    'Sb': 'Sb/POTCAR',
    'Te': 'Te/POTCAR',
    'I': 'I/POTCAR',
    'Xe': 'Xe/POTCAR',
    'Cs': 'Cs_sv/POTCAR',
    'Ba': 'Ba_sv/POTCAR',
    'La': 'La/POTCAR',
    'Ce': 'Ce_3/POTCAR',
    'Pr': 'Pr_3/POTCAR',
    'Nd': 'Nd_3/POTCAR',
    'Pm': 'Pm_3/POTCAR',
    'Sm': 'Sm_3/POTCAR',
    'Eu': 'Eu_3/POTCAR',
    'Gd': 'Gd_3/POTCAR',
    'Tb': 'Tb_3/POTCAR',
    'Dy': 'Dy_3/POTCAR',
    'Ho': 'Ho_3/POTCAR',
    'Er': 'Er_3/POTCAR',
    'Tm': 'Tm_3/POTCAR',
    'Yb': 'Yb_2/POTCAR',
    'Lu': 'Lu_3/POTCAR',
    'Hf': 'Hf_sv/POTCAR',
    'Ta': 'Ta_pv/POTCAR',
    'W': 'W_pv/POTCAR',
    'Re': 'Re_pv/POTCAR',
    'Os': 'Os_pv/POTCAR',
    'Ir': 'Ir/POTCAR',
    'Pt': 'Pt/POTCAR',
    'Au': 'Au/POTCAR',
    'Hg': 'Hg/POTCAR',
    'Tl': 'Tl_d/POTCAR',
    'Pb': 'Pb_d/POTCAR',
    'Bi': 'Bi/POTCAR',
    'Th': 'Th/POTCAR',
    'Pa': 'Pa/POTCAR',
    'U': 'U/POTCAR',
    'Np': 'Np/POTCAR',
    'Pu': 'Pu/POTCAR',
# 添加其他元素和对应的赝势路径
}

for element, relative_path in pseudo_paths.items():
    full_path = f"{base_pseudo_path}{relative_path}"
    pseudo_paths[element] = full_path


# 读取 POSCAR 文件
def read_poscar(poscar_path):  
    with open(poscar_path, 'r') as f:  
        lines = f.readlines()  
        elements = lines[5].split()
        return elements
# 创建大的 POTCAR 文件
def create_big_potcar(elements, pseudo_paths, output_path):
    with open(output_path, 'w') as f_out:
        for element in elements:
            if element in pseudo_paths:
                pseudo_path = pseudo_paths[element]
                with open(pseudo_path, 'r') as f_pseudo:
                    f_out.write(f_pseudo.read())
            else:
                print(f"赝势文件不存在或未定义：{element}")
def copypot(elments,pseudo_paths):
    import os 
    for element in elements:
        output_path2 = "POT-"+element
        if element in pseudo_paths:
            pseudo_path = pseudo_paths[element]
            with open(output_path2, 'w') as f_out:
                with open(pseudo_path, 'r') as f_pseudo:
                    f_out.write(f_pseudo.read())
        else:
            print(f"赝势文件不存在或未定义：{element}")
            
if __name__ == "__main__":
    poscar_path = "POSCAR"  # 输入文件名
    output_potcar_path = "POTCAR"  # 输出的大 POTCAR 文件名
    elements = read_poscar(poscar_path)
    create_big_potcar(elements, pseudo_paths, output_potcar_path)
    copypot(elements,pseudo_paths)
    print("大的 POTCAR 文件已创建")
