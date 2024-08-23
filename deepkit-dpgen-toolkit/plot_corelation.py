import os
import sys
from dpdata import System, LabeledSystem, MultiSystems
import matplotlib.pyplot as plt
import numpy as np

def prepare_all_compositions(testdata_path, pbfile_path):
    subdirectories = [os.path.join(testdata_path, d) for d in os.listdir(testdata_path) if os.path.isdir(os.path.join(testdata_path, d))]
    total_datas = MultiSystems()
    dft_energies = []
    ml_energies  = []
    dft_forces = []
    ml_forces  = []
    dft_virials = []
    ml_virials  = []
    for subdir in subdirectories:
        print(subdir)
        test_systems=LabeledSystem(subdir,  fmt="deepmd/npy")
        predict = test_systems.predict(pbfile_path)
        dft_energies.extend(test_systems['energies'])
        ml_energies.extend(predict['energies'])
        dft_forces.extend(test_systems['forces'])
        ml_forces.extend(predict['forces'])
        dft_virials.extend(test_systems['virials'])
        ml_virials.extend(predict['virials'])

    dft_forces   = np.concatenate(dft_forces, axis=0)
    ml_forces    = np.concatenate(ml_forces, axis=0)
    dft_virials  = np.concatenate(dft_virials, axis=0)
    ml_virials   = np.concatenate(ml_virials, axis=0)

    dft_forces   = dft_forces.flatten()
    ml_forces    = ml_forces.flatten()
    dft_virials  = dft_virials.flatten()
    ml_virials   = ml_virials.flatten()

    return dft_energies, ml_energies, dft_forces, ml_forces, dft_virials, ml_virials

def prepare_single_compositions(testdata_path, pbfile_path):
    test_systems = LabeledSystem(testdata_path, fmt = "deepmd/npy")
    predict = test_systems.predict(pbfile_path)
    dft_energies = test_systems['energies']
    ml_energies  = predict['energies']
    dft_forces = test_systems['forces']
    ml_forces  = predict['forces']
    dft_virials = test_systems['virials']
    ml_virials  = predict['virials']

    dft_forces   = dft_forces.flatten()
    ml_forces    = ml_forces.flatten()
    dft_virials  = dft_virials.flatten()
    ml_virials   = ml_virials.flatten()

    return dft_energies, ml_energies, dft_forces, ml_forces, dft_virials, ml_virials



if __name__ == "__main__":
    print("You can run it by: \npython plot_corelation.py ../../../../1.dp-data/2.testset/Ce1Sc2H22/ frozen_model.pb \nor \npython plot_corelation.py ../../../../1.dp-data/2.testset frozen_model.pb")
    # 制定训练集的路径
    testdata_path = os.path.abspath(sys.argv[1])
    # 制定模型路径
    pbfile_path = os.path.abspath(sys.argv[2])

    dft_energies, ml_energies = None, None
    if os.path.basename(testdata_path) == "2.testset":
        (dft_energies,
        ml_energies,
        dft_forces,
        ml_forces,
        dft_virials,
        ml_virials) = prepare_all_compositions(testdata_path, pbfile_path)
    else:
        (dft_energies,
        ml_energies,
        dft_forces,
        ml_forces,
        dft_virials,
        ml_virials) = prepare_single_compositions(testdata_path, pbfile_path)

    # 创建一个包含3个子图的Figure，每个子图占据一行中的一列  
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1行3列，调整figsize以适应你的需求  
    
    # 第一个子图：DFT能量 vs MPL预测能量  
    axs[0].scatter(dft_energies, ml_energies)  
    x_range = np.linspace(axs[0].get_xlim()[0], axs[0].get_xlim()[1])  
    axs[0].plot(x_range, x_range, 'r--', linewidth=0.25)  
    axs[0].set_xlabel("Energy of DFT")  
    axs[0].set_ylabel("Energy predicted by MPL")  
    axs[0].set_title("DFT vs MPL Energy")  # 可选：添加标题  
    
    # 第二个子图：DFT力 vs MPL预测力  
    axs[1].scatter(dft_forces, ml_forces)  # 注意：这里可能需要处理多维数据（如取模长或特定分量）  
    # 假设forces是三维的，这里我们仅作为示例取第一个分量的比较  
    # 如果forces是一维的，则直接使用  
    x_range = np.linspace(axs[1].get_xlim()[0], axs[1].get_xlim()[1])  
    axs[1].plot(x_range, x_range, 'r--', linewidth=0.25)  
    axs[1].set_xlabel("Force of DFT")  # 假设我们比较x分量  
    axs[1].set_ylabel("Force predicted by MPL")  
    axs[1].set_title("DFT vs MPL Force")  # 可选：添加标题  
    
    # 第三个子图：DFT维里系数 vs MPL预测维里系数  
    # 注意：维里系数的处理可能也需要根据你的数据维度进行调整  
    axs[2].scatter(dft_virials, ml_virials)  # 同样，这里可能需要处理多维数据  
    x_range = np.linspace(axs[2].get_xlim()[0], axs[2].get_xlim()[1])  
    axs[2].plot(x_range, x_range, 'r--', linewidth=0.25)  
    axs[2].set_xlabel("Virial of DFT")  
    axs[2].set_ylabel("Virial predicted by MPL")  
    axs[2].set_title("DFT vs MPL Virial")  # 可选：添加标题  
    
    # 显示整个Figure  
    plt.tight_layout()  # 自动调整子图参数, 使之填充整个图像区域  
    plt.show()  
    
    # 保存整个Figure到一个文件中  
    plt.savefig("DFT_vs_MPL_comparison.png")

