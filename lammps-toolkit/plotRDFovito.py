#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
#plt.tight_layout()


def smear(data, sigma):
    """
    Apply Gaussian smearing to spectrum y value.

    Args:
        sigma: Std dev for Gaussian smear function
    """ 
    # print(data[0, 0 + 1], data[0, 0]); print(np.shape(data)[0]); print(data); input()
    diff = [data[0, i + 1] - data[0, i] for i in range(np.shape(data)[0] - 1)]
    avg_x_per_step = np.sum(diff) / len(diff)
    # print(f"diff={diff}, avg_x_per_step={avg_x_per_step}, sigma / avg_x_per_step={sigma / avg_x_per_step}")
    data[1, :] = gaussian_filter1d(data[1, :], sigma / avg_x_per_step)
    return data

def plot_RDF(RDFdata, sigma, filename=None, tag=None):
    """ plot PXRD """
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
    datax = smear(RDFdata, sigma)
    #print("self.RDF",self.RDF)
    #plt.plot(datax[0, :], datax[1, :])
    plt.plot(RDFdata[0, :],RDFdata[1, :],label=tag,linewidth=2.0)
    plt.xlabel("$r(\AA)$",fontsize=12)
    plt.ylabel("$g(r)$",fontsize=12)
    if filename is None:
        plt.show()
        #pass
    else:
        plt.savefig(filename)
        plt.close()
    
def plot_PRDF(PRDFdata, sigma, filename=None, tag=None):
    """ plot PXRD """
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
    datax = smear(PRDFdata, sigma)
    #print("self.RDF",self.RDF)
    #plt.plot(datax[0, :], datax[1, :])
    for i in range(1, PRDFdata.shape[0]):
        print(df.columns[i])
        plt.plot(PRDFdata[0, :], PRDFdata[i, :], linewidth=2.0, label=df.columns[i])
    plt.xlabel("$r(\AA)$",fontsize=12)
    plt.ylabel("$g(r)$",fontsize=12)
    plt.legend()

    if filename is None:
        plt.show()
        #pass
    else:
        plt.savefig(filename)
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot RDF and PRDF")
    parser.add_argument('-i', '--input_file', type=str, default='rdf.txt', help="input just like RDF.txt ")
    parser.add_argument('-s', '--sigma', type=float, default=0.1, help="gauss broading")
    parser.add_argument("-p", "--plot", dest="plot", default=None, help="generate the plot to file, default: None")
    args = parser.parse_args()

    # 打开文件并读取前两行
    with open(args.input_file, 'r') as file:
        lines = [next(file) for _ in range(2)]

    # 第二行包含列名，去除注释符号 '#' 并分割成列名
    # lines[1] = "Pair  separation  distance" Cl-Cl Cl-Na Na-Na
    # lines[1].strip().lstrip('#').split()[3:] = ['Cl-Cl', 'Cl-Na', 'Na-Na']
    column_names = ['Pair  separation  distance']+ lines[1].strip().lstrip('#').split()[3:]

    # 读取剩余的数据，跳过前两行，并使用刚才提取的列名
    df = pd.read_table(args.input_file, comment='#', sep='\s+', skiprows=2, names=column_names)

    if df.shape[1] == 2:
        plot_RDF(df.to_numpy().T, args.sigma, args.plot)
    else:
        plot_PRDF(df.to_numpy().T, args.sigma, args.plot)
