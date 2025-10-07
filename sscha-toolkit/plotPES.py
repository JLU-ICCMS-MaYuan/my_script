import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# 设置随机种子以确保可重复性
np.random.seed(42)

# 创建网格
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)

# 1. 修改后的经典势能面 V(R) - 降低崎岖程度
def classical_potential(x, y):
    # 减少高斯峰的数量和强度
    z = (0.5*np.sin(0.8*x)*np.cos(0.8*y) + 
         0.3*np.exp(-((x-7)**2 + (y-7)**2)/3) -
         0.3*np.exp(-((x+2)**2 + (y+2)**2)/4) +
         0.2*np.exp(-((x-3)**2 + (y+3)**2)/6) -
         0.2*np.exp(-((x+3)**2 + (y-3)**2)/8))
    return z

# 2. 量子效应势能面 E(R) - 保持原有平滑度
def quantum_potential(x, y):
    # 单个全局最小值
    z = 0.5*(x**2 + y**2) * (1 + 0.1*np.cos(0.5*x)*np.sin(0.5*y))
    return z

# 计算势能面
Z_classical = classical_potential(X, Y)
Z_quantum = quantum_potential(X, Y)

# 创建图形
fig = plt.figure(figsize=(12, 10))

# 调整布局
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.3)

# 1. 绘制修改后的经典势能面
ax1 = fig.add_subplot(211, projection='3d')
surf1 = ax1.plot_surface(X, Y, Z_classical, cmap=cm.rainbow, 
                        linewidth=0, antialiased=True, alpha=0.8)
ax1.set_title('Classical description V(R) (Reduced Ruggedness)', color='brown', y=1.02, fontsize=12)
ax1.set_zlabel('Energy', fontsize=10)
ax1.view_init(elev=30, azim=-60)

# 2. 绘制量子效应势能面
ax2 = fig.add_subplot(212, projection='3d')
surf2 = ax2.plot_surface(X, Y, Z_quantum, cmap=cm.rainbow, 
                        linewidth=0, antialiased=True, alpha=0.8)
ax2.set_title('with Quantum effects E(R)', color='brown', y=1.02, fontsize=12)
ax2.set_zlabel('Energy', fontsize=10)
ax2.view_init(elev=30, azim=-60)

# 添加颜色条
fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10, pad=0.1)
fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10, pad=0.1)

# 调整视角和显示
plt.tight_layout()
plt.show()