#!/usr/bin/env python3
import argparse
import os
import re
import shutil

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar


def get_struct_info(struct_path):
    '''
    input parameter:
        struct:         pymatgen class structure
        output_poscar:  the output path of primitive cell
    return:
        None
    '''
    atom   = read(struct_path)
    struct = AseAtomsAdaptor.get_structure(atom)
    spa = SpacegroupAnalyzer(struct)
    # bstruct = spa.get_conventional_standard_structure()
    # pstruct = spa.get_primitive_standard_structure()
    # Poscar(pstruct).write_file(output_poscar.joinpath("PPOSCAR"))

    # 处理PPOSCAR的pymatgen对象
    # 获得元素名称 和 每种元素的原子个数
    composition        = struct.composition.get_el_amt_dict()
    species            = struct.composition.elements
    # 获得体系的 化学配比
    system_name        = struct.composition.formula.replace(" ", "")
    # 获得元素种类的个数
    species_quantity   = len(composition)
    # 获得每种元素的相对原子质量
    all_atoms_quantity = int(sum(composition.values()))

    return system_name

def get_q_from_dyn0(dyn0_path):
    '''
    input  : directory        dyn0文件的路径
             
    return : self.qtot 总q点数
             self.qirreduced 不可约q点数
             self.qirreduced_coords 不可约q点坐标
    '''
    if not os.path.exists(dyn0_path):
        raise FileExistsError ("dyn0 file doesn't exist!")
    content = open(dyn0_path, "r").readlines()
    # check qtot is right or not! 
    qpoints   = list(map(int, content[0].strip("\n").split()))
    qtot_list = list(map(int, qpoints))
    qtot      = qtot_list[0] * qtot_list[1] * qtot_list[2]

    def find_q(item):
        if re.search(r"E[\+|\-]", item):
            return item

    qirreduced = int(content[1])
    _q_coordinate_list = list(filter(find_q, content))
    qirreduced_coords = [q_string.strip("\n").split() for q_string in _q_coordinate_list]
    
    return qtot, qirreduced, qirreduced_coords

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="指定合并tmp文件时的方式")
    parser.add_argument("-s", "--specifyQ_type", type=str, required=True, help="只有by_index和by_coords")
    args = parser.parse_args()
    specifyQ_type   = args.specifyQ_type

    system_name     = get_struct_info("scffit.out")
    dyn0_path       = system_name+".dyn0"
    print(dyn0_path)
    qtot, qirreduced, qirreduced_coords = get_q_from_dyn0(dyn0_path)
    
    #准备好合并后的tmp
    merged_tmp_path = os.path.abspath("tmp")
    merged_ph0_path = os.path.join(merged_tmp_path, "_ph0")
    merged_phsave_path = os.path.join(merged_ph0_path, system_name+".phsave")
    if (not os.path.exists(merged_tmp_path)):
        os.makedirs(merged_tmp_path)
    if (not os.path.exists(merged_ph0_path)):
        os.makedirs(merged_ph0_path)
    if (not os.path.exists(merged_phsave_path)):
        os.makedirs(merged_phsave_path)
        
    if specifyQ_type == "by_index":    
        for i in range(qirreduced):
            print(i+1)
            if i==0 :
                splited_dv1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf1")
                merged_dv1_path     = os.path.join(merged_ph0_path,         system_name+".dvscf1")
                shutil.copy(splited_dv1_path, merged_dv1_path)

                splited_dv_paw1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf_paw1")
                merged_dv_paw1_path     = os.path.join(merged_ph0_path,         system_name+".dvscf_paw1")
                shutil.copy(splited_dv_paw1_path, merged_dv_paw1_path)
    
                splited_control_ph_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave", "control_ph.xml")
                merged_control_ph_path  = os.path.join(merged_ph0_path,         system_name+".phsave", "control_ph.xml")
                shutil.copy(splited_control_ph_path, merged_control_ph_path)
            else:
                q_num_path          = os.path.join(merged_ph0_path,         system_name+".q_"+str(i+1))
                if not os.path.exists(q_num_path):
                    os.makedirs(q_num_path)
                splited_dv1_path1    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".q_"+str(i+1), system_name+".dvscf1")
                merged_dv1_path      = os.path.join(q_num_path,                                          system_name+".dvscf1")    
                shutil.copy(splited_dv1_path, merged_dv1_path)
      
                splited_dv_paw1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".q_"+str(i+1), system_name+".dvscf_paw1")
                merged_dv_paw1_path     = os.path.join(q_num_path,                                          system_name+".dvscf_paw1")
                shutil.copy(splited_dv_paw1_path, merged_dv_paw1_path)
    
            splited_phsave_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave")
            merged_phsave_path  = os.path.join(merged_ph0_path, system_name+".phsave")
            for src in os.listdir(splited_phsave_path):
                print(src)
                src_path = os.path.join(splited_phsave_path, src)
                if "dynmat" in src_path or "elph" in src_path:
                    shutil.copy(src_path, merged_phsave_path)
    
            splited_patterns_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave", "patterns."+str(i+1)+".xml")
            merged_patterns_path  = os.path.join(merged_ph0_path,         system_name+".phsave", "patterns."+str(i+1)+".xml")
            shutil.copy(splited_patterns_path, merged_patterns_path)
    if specifyQ_type == "by_coords":    
        for i in range(qirreduced):
            print(i+1)
            if i==0 :
                splited_dv1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf1")
                merged_dv1_path     = os.path.join(merged_ph0_path,         system_name+".dvscf1")
                shutil.copy(splited_dv1_path, merged_dv1_path)
    
                splited_control_ph_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave", "control_ph.xml")
                merged_control_ph_path  = os.path.join(merged_ph0_path,         system_name+".phsave", "control_ph.xml")
                shutil.copy(splited_control_ph_path, merged_control_ph_path)

                splited_dv_paw1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf_paw1")
                merged_dv_paw1_path     = os.path.join(merged_ph0_path,         system_name+".dvscf_paw1")
                shutil.copy(splited_dv_paw1_path, merged_dv_paw1_path)
            else:
                q_num_path          = os.path.join(merged_ph0_path,         system_name+".q_"+str(i+1))
                if not os.path.exists(q_num_path):
                    os.makedirs(q_num_path)
                splited_dv1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf1")
                merged_dv1_path      = os.path.join(q_num_path,              system_name+".dvscf1")    
                shutil.copy(splited_dv1_path, merged_dv1_path)
      
                splited_dv_paw1_path    = os.path.join(str(i+1), "tmp", "_ph0", system_name+".dvscf_paw1")
                merged_dv_paw1_path     = os.path.join(merged_ph0_path,         system_name+".dvscf_paw1")
                shutil.copy(splited_dv_paw1_path, merged_dv_paw1_path)
    
            splited_phsave_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave")
            merged_phsave_path  = os.path.join(merged_ph0_path, system_name+".phsave")
            for src in os.listdir(splited_phsave_path):
                print(src)
                src_path = os.path.join(splited_phsave_path, src)
                if "dynmat" in src_path or "elph" in src_path:
                    shutil.copy(src_path, merged_phsave_path)
    
            splited_patterns_path = os.path.join(str(i+1), "tmp", "_ph0", system_name+".phsave", "patterns.1.xml")
            merged_patterns_path  = os.path.join(merged_ph0_path,         system_name+".phsave", "patterns."+str(i+1)+".xml")
            shutil.copy(splited_patterns_path, merged_patterns_path)
            os.system(f"sed -i '4s/1/{str(i+1)}/' {merged_patterns_path}")
