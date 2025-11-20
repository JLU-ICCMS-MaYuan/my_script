#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# 设置图像全局属性
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

# 读取并绘制每一张图的函数
def plot_subplot(ax, filename, xlabel, ylabel, logscale=True):
    data = np.loadtxt(filename)
    x, y = data[:, 0], data[:, 1]
    if logscale:
        ax.set_yscale('log')
    ax.plot(x, y, 'o', markersize=5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls="--", lw=0.5)

# 创建 2x2 子图
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

plot_subplot(axs[0, 0], "decay.H", r"|R_e| (\u00C5)", r"max $H_{nm}$ (Ry)")
plot_subplot(axs[0, 1], "decay.dynmat", r"|R_p| (\u00C5)", r"max $D_{nm}$ (Ry)")
plot_subplot(axs[1, 0], "decay.epmate", r"|R_e| (\u00C5)", r"max $g_{nm}$ (Ry)")
plot_subplot(axs[1, 1], "decay.epmatp", r"|R_p| (\u00C5)", r"max $g_{nm}$ (Ry)")

plt.tight_layout()
plt.savefig("decay.png", dpi=300)
plt.close()

