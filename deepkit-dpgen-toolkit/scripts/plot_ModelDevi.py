import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def parepare_single_model_deviation(file_path):
    # 读取文件
    steps = []
    max_devi_f = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        # 跳过以 # 开头的注释行
        if line.startswith('#'):
            continue
        # 分割行内容，按空格分割
        #print(line)
        parts = line.split()
        if parts:
            # print(parts)
            # 解析数据，添加到列表
            steps.append(int(parts[0]))
            max_devi_f.append(float(parts[4]))
    return max_devi_f

def prapare_accuracies_from_dpgenlog():
    content = os.popen('grep -B1 -P "iter\.\d+ task 07" dpgen.log | grep "summary  accurate_ratio"').readlines()
    pattern = re.compile(r'summary\s+accurate_ratio:\s+([\d.]+)%\s+candidata_ratio:\s+([\d.]+)%\s+failed_ratio:\s+([\d.]+)%\s+in\s+(\d+)\s+structures')
    accurate_ratios = []; candidata_ratios = []; failed_ratios = []; structuress = []
    for cont in content:
        accurate_ratio, candidata_ratio, failed_ratio, structures = pattern.search(cont).groups()
        print(f"Accurate Ratio: {accurate_ratio}%   Candidata Ratio: {candidata_ratio}%   Failed Ratio: {failed_ratio}%   Structures: {structures}")
        accurate_ratios.append(float(accurate_ratio))
        candidata_ratios.append(float(candidata_ratio))
        failed_ratios.append(float(failed_ratio))
        structuress.append(structures)
        
    return accurate_ratios, candidata_ratios, failed_ratios, structuress
   
    

def plot_single_model_deviation(max_devi_f: list[float], ax=None):

    filtered_max_devi_f = [x for x in max_devi_f if 0 <= x <= 5]  
    if ax is None:
        fig, ax = plt.subplots()
    ax.hist(filtered_max_devi_f, bins=100, color='blue', edgecolor='black', alpha=0.7)
    ax.set_title('Distribution of Max Devi F')
    ax.set_xlabel('Max Devi F')
    ax.set_ylabel('Frequency')
    ax.grid(True)

    plt.savefig("max_devi_f_distribution.png")
    plt.show()


def plot_many_model_deviations(max_devi_fs:list[list[float]], accuracies: list[float]):
    num_models = len(max_devi_fs)
    print(num_models)
    fig, axes = plt.subplots(1, num_models, figsize=(4 * num_models, 8), sharey=True) #  sharey=True 共用y轴
    
    for i, max_devi_f in enumerate(max_devi_fs):
        filtered_max_devi_f = [x for x in max_devi_f if 0 <= x <= 5]  
        counts, bins = np.histogram(filtered_max_devi_f, bins=100)
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
        axes[i].axhline(y=2, color='red', linestyle='--')
        axes[i].axhline(y=1, color='red', linestyle='--')
        axes[i].set_title(f'Iteration No. {i+1}')

        # 设置横坐标的间隔为100
        max_x_value = max(counts)
        axes[i].set_xticks(range(0, max_x_value + 100, 200))

        # 只保留第一个axes的y轴标签
        if i > 0:
            axes[i].set_ylabel('')

    # 获取第一个子图的高度和底部位置
    pos0 = axes[0].get_position()  # 返回子图的 [left, bottom, width, height]

    # 动态计算新轴的位置和大小
    accuracy_ax_left_margin  = axes[0].get_position().x0  # 第0个图的左边位置
    accuracy_ax_right_margin = axes[-1].get_position().x1 # 最后一个图的右边位置
    accuracy_ax_bottom_margin= axes[0].get_position().y0  # 图的底部位置
    accuracy_ax_width = accuracy_ax_right_margin - accuracy_ax_left_margin
    accuracy_ax_height = axes[0].get_position().y1 - axes[0].get_position().y0 #  y0 图的底部位置 y1图的顶部位置
    accuracy_ax = fig.add_axes([accuracy_ax_left_margin, accuracy_ax_bottom_margin, accuracy_ax_width, accuracy_ax_height], frameon=False) # 新的轴覆盖整个子图区域
    
    accuracy_ax.plot(range(1, num_models + 1), accuracies, 'go--', label='Accuracy')
    accuracy_ax.legend()
    accuracy_ax.set_xlim(0.5, num_models + 0.5)
    accuracy_ax.set_ylim(0, 100)  # 假设精确度范围为 0 到 1


    accuracy_ax.set_ylabel('Accurate Rate (%)')
    accuracy_ax.yaxis.set_label_position('right')

    accuracy_ax.set_xlabel('Frequency')
    accuracy_ax.xaxis.set_label_position('bottom')

    accuracy_ax.xaxis.set_ticks_position('top')  # 将 x 轴刻度放到顶部
    accuracy_ax.xaxis.set_label_position('top')  # 将 x 轴标签放到顶部
    
    accuracy_ax.yaxis.tick_right()        # 将 y 轴的刻度标签放到右侧
    accuracy_ax.xaxis.set_visible(False)  # 隐藏 x 轴，因为它会和子图的 x 轴重叠

    

    # 调整图的水平距离为0
    plt.subplots_adjust(wspace=0.1) 
    #plt.tight_layout()
    plt.savefig("max_devi_fs_distribution.png")
    plt.show()


    # 调整图的水平距离为0
    plt.subplots_adjust(wspace=0) 
    #plt.tight_layout()
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
        pattern = re.compile(r'^iter\.\d+$')
        # 筛选出符合条件的文件或文件夹
        matching_files = [p for p in os.listdir() if pattern.match(p)]
        #print(matching_files)
        # 根据 'iter.' 后面的数字进行排序
        matching_files.sort(key=lambda x: int(x.split('.')[-1]))
        #print(matching_files)
        for p in os.listdir():
            model_devi_path = os.path.join(file_path, p, "01.model_devi/model_devi_results/Model_Devi.out")
            if os.path.exists(model_devi_path): 
                print(model_devi_path, '----------------------')
                max_devi_f = parepare_single_model_deviation(model_devi_path)
                max_devi_fs.append(max_devi_f)
        accurate_ratios, candidata_ratios, failed_ratios, structuress = prapare_accuracies_from_dpgenlog()
        print(accurate_ratios)
        plot_many_model_deviations(max_devi_fs, accurate_ratios)