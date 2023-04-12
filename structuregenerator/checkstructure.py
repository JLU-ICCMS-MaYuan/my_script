from pymatgen.core.structure import Structure

from structuregenerator.clathrate import HydrideClathrate

def checkclathrate(pmg_struct: Structure):
    clathrate = HydrideClathrate(pmg_struct)
    # if clathrate.remain_H_ratio > 0.44:
    #     return False
    # if clathrate.fraction_of_hydrogen_volume < 0.3:
    #     return False
    # if clathrate.shr_num_avg < 1.6:
    #     return False
    if clathrate.cage_regularity_avg > 0.05:
        return False
    if clathrate.h2h_network_regularity_avg > 0.15:
        return False
    if clathrate.h2h_1nearest[0][0] < 0.9:
        return False
    if clathrate.nonh2nonh_1nearest[0][0] < 2.2:
        return False
    
    return True


if __name__ == "__main__":
    import sys
    from pathlib import Path

    print("如果你想使用这个脚本, 你需要指定一个路径, 这个路径下存放着所有待判断的结构")
    clathrates = []
    for file in Path("/work/home/may/workplace/6.new_hydride/1.test-my-method/2.lah10-my-method/1.calypso-structure/100").glob("CONTCAR_*"):
        filename = file.name
        print(filename)
        struct = Structure.from_file(file)
        if checkclathrate(struct):
            clathrates.append(filename)
    
    print("最终结果是：")
    for name in clathrates:
        print(name)
    print("总共数量是：", len(clathrates))