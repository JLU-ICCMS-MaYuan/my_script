import matplotlib.pyplot as plt
import numpy as np
import sys,os

def plot_files(file_names):
    # 读取输入文件的数据
    data = []
    for file_name in file_names:
        data.append(np.loadtxt(file_name))

    # 绘制曲线图
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    for i in range(len(data)):
        for j in range(1, data[i].shape[1]):
            if j == 1:
                plt.plot(data[i][:,0], data[i][:,j], color=colors[i], label=os.path.basename(file_names[i]))
            else:
                plt.plot(data[i][:,0], data[i][:,j], color=colors[i])

    # 添加图例
    plt.legend(loc='upper left')
    # 添加x轴和y轴标签
    plt.xlabel('snna_v3')
    plt.ylabel('Frequency(cm-1)')


    # 设置x轴范围
    plt.xlim(data[0][0,0], data[0][-1,0])

    # 保存图像
    plt.savefig('snna_v3.jpg')

if __name__ == '__main__':
    file_names = sys.argv[1:]
    plot_files(file_names)


