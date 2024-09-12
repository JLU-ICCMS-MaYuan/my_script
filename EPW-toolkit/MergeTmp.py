#!/usr/bin/env python3
import sys
import os
import re


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
    dyn0_path = sys.argv[1] 
    qtot, qirreduced, qirreduced_coords = get_q_from_dyn0(dyn0_path)

    #准备好合并后的tmp
    merged_tmp_path = os.path.abspath("tmp")
    merged_ph0_path = os.path.join(merged_tmp_path, "_ph0")
    merged_phsave_path = os.path.join(merged_tmp_path, "Nb4H14.phsave")
    if (not os.path.exists(merged_tmp_path)) and (not os.path.exists(merged_ph0_path)) and (not os.path.exists(merged_phsave_path)):
        os.makedirs(merged_tmp_path)
    for i in range(qirreduced):
        if i==0 :
            splited_ph0_path = os.path.join(str(i+1), "tmp", "Nb4H14.Nb4H14.dv1")
            splited_phsave_path = os.path.join(str(i+1), "tmp", "Nb4H14.phsave")
