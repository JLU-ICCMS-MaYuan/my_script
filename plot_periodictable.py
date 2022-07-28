#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
import matplotlib as mpl  # 导入Matplotlib库
import matplotlib.cm as cm  # 导入Matplotlib颜色配置模块
import matplotlib.patches as patches  # 导入Matplotlib形状绘制模块
import matplotlib.pyplot as plt  # 导入Matplotlib绘图模块
import mendeleev  # 导入元素周期表库（包含118种元素的基本性质）

# 绘制热力图数据
plot_data = {'O': 9, 'Te': 9, 'F': 9, 'S': 9, 'Na': 9, 'K': 9, 'N': 8, 'Li': 8,
             'I': 8, 'Rb': 8, 'Si': 7, 'Cd': 7, 'Cl': 7, 'Zn': 7, 'H': 7,
             'Bi': 7, 'Br': 7, 'P': 6, 'Sn': 6, 'Ca': 6, 'Au': 6, 'Al': 5,
             'As': 5, 'Ga': 5, 'C': 5, 'Ge': 5, 'Sr': 5, 'Se': 5, 'Be': 5,
             'B': 5, 'Cs': 5, 'Mg': 5, 'Ag': 5, 'Pb': 4, 'In': 4, 'Ti': 4,
             'Cu': 4, 'Zr': 4, 'Sb': 4, 'Tl': 4, 'Sc': 4, 'Y': 4, 'Hg': 4,
             'Ba': 4, 'La': 4, 'Hf': 4, 'Og': 1}

# 元素周期表中cell的设置
# cell的大小
cell_length = 1
# 各个cell的间隔
cell_gap = 0.1
# cell边框的粗细
cell_edge_width = 0.5

# 获取各个元素的原子序数、周期数（行数）、族数（列数）以及绘制数据（没有的设置为0）
elements = []
for i in range(1, 119):
    ele = mendeleev.element(i)
    ele_group, ele_period = ele.group_id, ele.period

    # 将La系元素设置到第8行
    if 57 <= i <= 71:
        ele_group = i - 57 + 3
        ele_period = 8
    # 将Ac系元素设置到第9行
    if 89 <= i <= 103:
        ele_group = i - 89 + 3
        ele_period = 9

    elements.append([i, ele.symbol, ele_group, ele_period,
                     plot_data.setdefault(ele.symbol, 0)])

# 设置La和Ac系的注解标签
elements.append([None, 'LA', 3, 6, None])
elements.append([None, 'AC', 3, 7, None])
elements.append([None, 'LA', 2, 8, None])
elements.append([None, 'AC', 2, 9, None])

# 新建Matplotlib绘图窗口
fig = plt.figure(figsize=(10, 5))
# x、y轴的范围
xy_length = (20, 11)

# 获取YlOrRd颜色条
my_cmap = cm.get_cmap('YlOrRd')
# 将plot_data数据映射为颜色，根据实际情况调整
norm = mpl.colors.Normalize(1, 10)
# 设置超出颜色条下界限的颜色（None为不设置，即白色）
my_cmap.set_under('None')
# 关联颜色条和映射
cmmapable = cm.ScalarMappable(norm, my_cmap)
# 绘制颜色条
plt.colorbar(cmmapable, drawedges=False)

# 绘制元素周期表的cell，并填充属性和颜色
for e in elements:
    ele_number, ele_symbol, ele_group, ele_period, ele_count = e
    print(ele_number, ele_symbol, ele_group, ele_period, ele_count)

    if ele_group is None:
        continue

    # x, y定位cell的位置
    x = (cell_length + cell_gap) * (ele_group - 1)
    y = xy_length[1] - ((cell_length + cell_gap) * ele_period)

    # 增加 La, Ac 系元素距离元素周期表的距离
    if ele_period >= 8:
        y -= cell_length * 0.5

    # cell中原子序数部位None时绘制cell边框并填充热力颜色
    # 即不绘制La、Ac系注解标签地边框以及颜色填充
    if ele_number:
        fill_color = my_cmap(norm(ele_count))
        rect = patches.Rectangle(xy=(x, y),
                                 width=cell_length, height=cell_length,
                                 linewidth=cell_edge_width,
                                 edgecolor='k',
                                 facecolor=fill_color)
        plt.gca().add_patch(rect)

    # 在cell中添加原子序数属性
    plt.text(x + 0.04, y + 0.8,
             ele_number,
             va='center', ha='left',
             fontdict={'size': 6, 'color': 'black', 'family': 'Helvetica'})
    # 在cell中添加元素符号
    plt.text(x + 0.5, y + 0.5,
             ele_symbol,
             va='center', ha='center',
             fontdict={'size': 9, 'color': 'black', 'family': 'Helvetica', 'weight': 'bold'})
    # 在cell中添加热力值
    plt.text(x + 0.5, y + 0.12,
             ele_count,
             va='center', ha='center',
             fontdict={'size': 6, 'color': 'black', 'family': 'Helvetica'})

# x, y 轴设置等比例（1:1）（使cell看起来是正方形）
plt.axis('equal')
# 关闭坐标轴
plt.axis('off')
# 裁剪空白边缘
plt.tight_layout()
# 设置x, y轴的范围
plt.ylim(0, xy_length[1])
plt.xlim(0, xy_length[0])

# 将图保存为*.svg矢量格式
plt.savefig('./periodic_table.svg')
# 显示绘图窗口
plt.show()