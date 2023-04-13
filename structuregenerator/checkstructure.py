from pymatgen.core.structure import Structure

from structuregenerator.clathrate import HydrideClathrate

def checkclathrate(
        pmg_struct: Structure,
        remain_H_ratio_UPPERSTD: float,
        fraction_of_hydrogen_volume_LOWERSTD: float,
        shr_num_avg_LOWERSTD: float,
        cage_regularity_avg_UPPERSTD: float,
        h2h_network_regularity_avg_UPPERSTD: float,
        h2h_1nearest_LOWERSTD: float,
        nonh2nonh_1nearest_LOWERSTD: float,
        ):
    clathrate = HydrideClathrate(pmg_struct)

    # 推荐设置
    # remain_H_ratio_UPPERSTD = 0.44
    # fraction_of_hydrogen_volume_LOWERSTD = 0.3
    # shr_num_avg_LOWERSTD = 1.6
    # cage_regularity_avg_UPPERSTD = 0.05
    # h2h_network_regularity_avg_UPPERSTD = 0.15
    # h2h_1nearest_LOWERSTD = 0.9
    # nonh2nonh_1nearest_LOWERSTD = 2.2
    if clathrate.remain_H_ratio > remain_H_ratio_UPPERSTD:
        return False
    if clathrate.fraction_of_hydrogen_volume < fraction_of_hydrogen_volume_LOWERSTD:
        return False
    if clathrate.shr_num_avg < shr_num_avg_LOWERSTD:
        return False
    if clathrate.cage_regularity_avg > cage_regularity_avg_UPPERSTD:
        return False
    if clathrate.h2h_network_regularity_avg > h2h_network_regularity_avg_UPPERSTD:
        return False
    if clathrate.h2h_1nearest[0][0] < h2h_1nearest_LOWERSTD:
        return False
    if clathrate.nonh2nonh_1nearest[0][0] < nonh2nonh_1nearest_LOWERSTD:
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