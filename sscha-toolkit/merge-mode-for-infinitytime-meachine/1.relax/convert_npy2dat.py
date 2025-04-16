import os
import argparse
import cellconstructor as CC
from sscha.Ensemble import Ensemble

def main(NQIRR, DYN_POP_IDX, END_POP_IDX):


    # 加载动力学数据
    dyn = CC.Phonons.Phonons(f"dyn_pop{DYN_POP_IDX}_", nqirr=NQIRR)
    ens = Ensemble(dyn0=dyn, T0=0, supercell=dyn.GetSupercell())

    # 加载并保存系综数据
    ens.load_bin(data_dir="ensembles", population_id=END_POP_IDX)
    ens.save("data_ensemble_manual", population=END_POP_IDX)

if __name__ == "__main__":
    # 设置 argparse 参数解析
    parser = argparse.ArgumentParser(description="处理声子动力学和系综数据")
    parser.add_argument("-n", "--nqirr", type=int, required=True, help="不可约 q 点的数量")
    parser.add_argument("-i", "--dyn_pop_idx", type=int, required=True, help="要读取的 dyn 文件的代编号")
    
    args = parser.parse_args()
    
    NQIRR = args.nqirr
    DYN_POP_IDX = args.dyn_pop_idx
    END_POP_IDX = DYN_POP_IDX + 1
    main(NQIRR, DYN_POP_IDX, END_POP_IDX)