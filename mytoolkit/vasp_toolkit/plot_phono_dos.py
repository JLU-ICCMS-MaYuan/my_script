#! /usr/bin/env python

import sys

import pandas as pd
import matplotlib.pyplot as plt


def load_figure(
    _max_frequency,
    _min_frequency,
    _max_qvector,
    _min_qvector,
    ):
    # 设置画布的大小6英寸*4英寸 与 清晰度 300dpi
    fig = plt.figure(figsize=(6,4), dpi=300)
    # 分割画布为1行2列，左右画布长度比为：4：1.5
    # 左右画布间隙大小为0.08
    gs = fig.add_gridspec(
        1, 2,
        width_ratios=(4,1.2),
        left=0.1, right=0.9, bottom=0.1, top=0.9,
        wspace=0.08, hspace=0.05,
    )
    # 创建spectrum的坐标系
    ax_phono = fig.add_subplot(gs[0, 0])
    # 处理phono坐标系中的x轴:        axis="x"
    # 不显示phono坐标系中的左侧刻度值  labelbottom=False
    # 不显示phono坐标系中的左侧刻度线  bottom=False
    ax_phono.tick_params(axis="x", bottom=False, labelbottom=False)

    # 处理phono坐标系中的y轴:            axis="y"
    # 显示phono坐标系中的y轴左侧刻度线朝里: direction="in"
    ax_phono.tick_params(axis="y", direction="in")


    # 创建dos的坐标系, 并于spectrum共享y轴坐标
    ax_dos = fig.add_subplot(gs[0, 1], sharey=ax_phono)

    # 处理dos坐标系中的x轴:            axis="x"
    # 显示dos坐标系中的x轴左侧刻度线朝里: direction="in"
    ax_dos.tick_params(axis="x", direction="in")

    # 处理dos坐标系中的y轴:        axis="y"
    # 不显示dos坐标系中的左侧刻度值  labelleft=False
    # 不显示dos坐标系中的左侧刻度线  left=False
    ax_dos.tick_params(axis="y", labelleft=False, left=False)

    ax_phono.set_title("phono_spectrum")
    ax_dos.set_title("phono_dos")

    # 设置phono纵坐标范围
    ax_phono.set_ylim(_min_frequency, _max_frequency+5)

    # 设置phono横坐标范围
    ax_phono.set_xlim(_min_qvector, _max_qvector)

    return ax_phono, ax_dos


def plot_phono_dos(
        ax_phono,
        ax_dos,
        _frequency,
        _qvector,
):
    zero_line = [0] * len(_qvector)
    ax_phono.plot(
        _qvector, zero_line,
        color='b',
        linewidth=0.4,
        linestyle='dashed',
    )

    # 开始画phono_spectrum图
    ax_phono.plot(
        _qvector, _frequency,
        color='r',
        linewidth=0.5,
        linestyle='solid',
    )


def plot_high_symmtry_points(
        ax_phono,
        high_symmetry_points,
        max_frequency,
        min_frequency,
):
    # 打开phono_spectrum图片的x坐标的标签 labelbottom=True
    # 关闭phono_spectrum图片的x坐标的刻度 bottom=False
    ax_phono.tick_params(axis="x", bottom=False, labelbottom=True, labelsize=10.0, grid_color='b', grid_linestyle='dashed')
    # 自动增加x轴的网格线
    ax_phono.grid(axis='x')
    # 在指定的坐标下high_symmetry_points, 设置对应的高对称点字母
    namelabel = ['G', 'X', 'M', 'G', 'R', 'X']
    ax_phono.set_xticks(high_symmetry_points)
    ax_phono.set_xticklabels(namelabel)

    # 手动绘制网格线
    # for qpoints in high_symmetry_points:
    #     qvectors = []; frequency = []
    #     qvectors.extend([qpoints, qpoints])
    #     frequency.extend([min_frequency, max_frequency+5])
    #
    #     ax_phono.plot(
    #         qvectors, frequency,
    #         color='b',
    #         linewidth=0.5,
    #         linestyle='dashed',
    #     )


if __name__ == "__main__":

    banddat_path = sys.argv[1]

    # 读入band.dat
    # index_col=0  以第0列作为index
    # header=None 不设置表头
    # sep='\s+'   用空格作为分隔符
    # skiprows=2  忽略前2行不读入
    df = pd.read_table(banddat_path, index_col=0, header=None, sep='\s+', skiprows=2)
    # 获取频率列表，最大频率，最小频率
    frequency = df.values
    max_frequency = df[1].max()
    min_frequency = df[1].min()
    # 获取经过投影后的一维q点坐标列表，获取最大坐标，获取最小坐标
    qvector = df.index.values
    max_qvector = df.index.max()
    min_qvector = df.index.min()
    # 加载画布
    ax_phono, ax_dos = load_figure(
        max_frequency, min_frequency,
        max_qvector, min_qvector,
    )

    # 开始绘制每一条声子谱
    cont = True
    qvectors = []
    frequencies = []
    with open(banddat_path, 'r') as f:
        title = f.readline()
        high_symmetry_points = f.readline()
        high_symmetry_points = list(map(float, high_symmetry_points.strip("\n").strip("#").split()))
        plot_high_symmtry_points(
            ax_phono,
            high_symmetry_points,
            max_frequency,
            min_frequency,
        )

        while cont:
            cont = f.readline()
            if cont:
                if cont != '\n':
                    data = cont.strip("\n").split()
                    qvectors.append(float(data[0]))
                    frequencies.append(float(data[1]))
                else:
                    if len(qvectors) == len(frequencies) == 51:
                        max_qvector = max(qvectors)
                        min_qvector = min(qvectors)
                        max_frequency = max(frequencies)
                        min_frequency = min(frequencies)
                        plot_phono_dos(
                            ax_phono, ax_dos,
                            frequencies, qvectors,
                        )
                        # 置空列表 qvector_frequency
                        qvectors = []
                        frequencies = []
                    # else:
                    #     print("cont = {}".format(cont))
                    #     print("now there are {} points!, it's not equal to 51".format(len(qvectors)))
            else:
                print("walk to the bottom of the file")

    # 保存图片
    plt.savefig("phono-dos.png") 