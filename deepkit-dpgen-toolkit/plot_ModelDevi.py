import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def parepare_single_model_deviation(file_path):
    # 读取文件
    steps = []
    max_devi_f = []
    with open(file_path, 'r') as file:
        for line in file:
            # 跳过以 # 开头的注释行
            if line.startswith('#'):
                continue
            # 分割行内容，按空格分割
            parts = line.split()
            # 解析数据，添加到列表
            steps.append(int(parts[0]))
            max_devi_f.append(float(parts[4]))
    return max_devi_f

def plot_single_model_deviation(max_devi_f: list[float], ax=None):

    if ax is None:
        fig, ax = plt.subplots()
    ax.hist(max_devi_f, bins=100, color='blue', edgecolor='black', alpha=0.7)
    ax.set_title('Distribution of Max Devi F')
    ax.set_xlabel('Max Devi F')
    ax.set_ylabel('Frequency')
    ax.grid(True)

    plt.savefig("max_devi_f_distribution.png")
    plt.show()

def plot_many_model_deviations(max_devi_fs:list[list[float]]):
    num_models = len(max_devi_fs)
    fig, axes = plt.subplots(1, num_models, figsize=(4 * num_models, 8), sharey=True) #  sharey=True 共用y轴
    plt.subplots_adjust(wspace=-0.5) 
    
    for i, max_devi_f in enumerate(max_devi_fs):
        counts, bins = np.histogram(max_devi_f, bins=100)
        # bins = 100 指定了直方图要分成多少个区间（bins）。在这里，你将数据分成了100个区间。如果 bins 是 [0, 1, 2, 3, 4]，那么这些区间是 [0, 1), [1, 2), [2, 3), [3, 4]。
        # counts 对应区间内的数据点个数
        print(len(counts), len(bins))
        axes[i].barh(bins[:-1], counts, height=bins[1] - bins[0], color='blue', edgecolor='black', alpha=0.7)
        # bins[:-1] 用于确定每个柱子在图中的纵坐标上的位置
        # height 别看叫height，其实是设置水平柱子的宽度。bins[1] - bins[0]的设置方法保证了水平放置的柱子之间在垂直方向上没有间隙
        axes[i].set_title('Distribution of Max Devi F')
        axes[i].set_xlabel('Frequency')
        axes[i].set_ylabel('Max Devi F')
        axes[i].grid(False)
        # 在 y=0.9 和 y=0.5 添加虚线
        axes[i].axhline(y=0.9, color='red', linestyle='--')
        axes[i].axhline(y=0.5, color='red', linestyle='--')
        axes[i].set_title(f'Model {i+1} Distribution')

        # 只保留第一个axes的y轴标签
        if i > 0:
            axes[i].set_ylabel('')

    plt.tight_layout()
    plt.savefig("max_devi_fs_distribution.png")
    plt.show()

if __name__ == "__main__":
    # 文件路径
    print('You can use it in three ways:')
    print('    1. Run it in `2.getdp-mod/iter.00000*/01.model_devi/model_devi_results`')
    print('    2. Run it in `2.getdp-mod` and do not need to specify any path')
    print('    3. Run it in anywhere and need to specify path of `2.getdp-mod`')
    if os.path.exists('Model_Devi.out'):
        file_path = 'Model_Devi.out'
        max_devi_f = parepare_single_model_deviation(file_path)
        plot_single_model_deviation(max_devi_f)
    else:
        if len(sys.argv) > 1:
            file_path = os.path.abspath(sys.argv[1])
        else:
            file_path = os.path.abspath(os.getcwd())
        max_devi_fs = []
        for p in os.listdir(file_path):
            model_devi_path = os.path.join(file_path, p, "01.model_devi/model_devi_results/Model_Devi.out")
            if os.path.exists(model_devi_path):
                max_devi_f = parepare_single_model_deviation(model_devi_path)
                max_devi_fs.append(max_devi_f)
        plot_many_model_deviations(max_devi_fs)
