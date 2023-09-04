#!/usr/bin/env python3

# 学习参考连接
    # https://mp.weixin.qq.com/s/y7if55laqS-R16It-QGUiw
    # 
import sys
import numpy as np

import matplotlib.pyplot as plt
# from matplotlib.collections import LineCollection
# from matplotlib.gridspec import GridSpec


from pymatgen.io.vasp.outputs import Vasprun
# from pymatgen.electronic_structure.core import Spin, OrbitalType
# from pymatgen.electronic_structure import plotter
from pymatgen.electronic_structure.plotter import BSDOSPlotter, BSPlotter, BSPlotterProjected, DosPlotter

eledosfile = sys.argv[1] # such as, eledos/vasprun.xml
ebandfile = sys.argv[2] # such as, eband/vasprun.xml
ebandkpoints = sys.argv[3] # such as, eband/KPOINTS

dosrun = Vasprun(eledosfile, parse_potcar_file=False)
dos_data = dosrun.complete_dos

# bands
bs_run = Vasprun(ebandfile, parse_projected_eigen=True, parse_potcar_file=False)
bs_data = bs_run.get_band_structure(
                                kpoints_filename=ebandkpoints,
                                # efermi=dosrun.efermi,
                                line_mode=True,
                                )


# 能带图加上态密度
# 原始的能带和态密度
plot1 = BSDOSPlotter(bs_projection = None, dos_projection= None)
plot1.get_plot(bs=bs_data,dos=dos_data)
plt.savefig('1.band_dos.png',dpi=300)
plt.show() 

# 将能带和态密度投影到元素
plot2 = BSDOSPlotter(bs_projection = 'elements', dos_projection= 'elements')
plot2.get_plot(bs=bs_data,dos=dos_data)
plt.savefig('2.band投影元素_dos投影元素.png',dpi=300)
plt.show() 

# 将能带投影到元素，态密度投影到轨道
plot3 = BSDOSPlotter(bs_projection = 'elements', dos_projection='orbitals')
plot3.get_plot(bs=bs_data,dos=dos_data)
plt.savefig('3.band投影元素_dos投影轨道.png',dpi=300)
plt.show() 

# 只画能带图
# 原始的能带
bsplot=BSPlotter(bs=bs_data)
bsplot.get_plot()
bsplot.save_plot("4.band.png", img_format='png')

# 将能带投影到元素
bsplotproj=BSPlotterProjected(bs=bs_data)
plt=bsplotproj.get_elt_projected_plots()
#plt = bsplotproj.get_plt_projected_plots_color()
# bsplotproj.save_plot("5.band投影元素.png")
plt.savefig("5.band投影元素.png")

# 将能带投影到不同元素的不同轨道，即角量子数
bsplotproj=BSPlotterProjected(bs=bs_data)
plt=bsplotproj.get_projected_plots_dots(
    dictio={'La':['s', 'p', 'd'],
            # 'Ce' :['s', 'p', 'd'],
            'Be':['s', 'p'],
            'H' :['s'],})
# bsplotproj.save_plot("6.band投影指定元素的指定角量子轨道.png")
plt.savefig("6.band投影指定元素的指定角量子轨道.png")

# 将能带投影到不同元素的不同轨道的分量，即磁量子数
bsplotproj=BSPlotterProjected(bs=bs_data)
plt=bsplotproj.get_projected_plots_dots_patom_pmorb(
    dictio={'Be':['px', 'py', 'pz'],
            'La':['px', 'py', 'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']},
    dictpa={'Be':[3, 4],
            'La':[1]})
# bsplotproj.save_plot("7.band投影指定元素的指定磁量子轨道和指定原子轨道.png")
plt.savefig("7.band投影指定元素的指定磁量子轨道和指定原子轨道.png")

# sum_atoms可以画同种原子的总投影，sum_morbs可以画同种原子不同轨道的总投影
plt=bsplotproj.get_projected_plots_dots_patom_pmorb(
    dictio={'La':['px', 'py', 'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2'],
            'Be':['px', 'py', 'pz']},
    dictpa={'Be':[3, 4],
            'La':[1]},
    sum_morbs={'La':['px', 'py', 'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2'],
               'Be':['px', 'py', 'pz']})
# bsplotproj.save_plot("8.band投影同种原子的轨道.png")
plt.savefig("8.band投影同种原子的轨道.png")



# 只画态密度
# 总态密度: stack表示是否填充颜色，sigma表示是否展宽
dostotal=DosPlotter(stack=False, sigma=0.5)
dostotal.add_dos('total dos', dos=dos_data)
dostotal.save_plot('9.dos总态密度.png', img_format=u'png')


# 投影态密度到轨道
dostotal=DosPlotter(stack=False, sigma=0.5)
dostotal.add_dos('total dos', dos=dos_data)
dostotal.add_dos_dict(dos_data.get_spd_dos())
dostotal.save_plot('10.dos投影态密度到轨道.png', img_format=u'png')

# 投影态密度到元素
dostotal=DosPlotter(stack=False, sigma=0.5)
dostotal.add_dos('total dos', dos=dos_data)
dostotal.add_dos_dict(dos_data.get_element_dos())
dostotal.save_plot('11.dos投影态密度到元素.png', img_format=u'png')