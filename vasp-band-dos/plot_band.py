from pymatgen.io.vasp import Vasprun, Poscar
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
import matplotlib.pyplot as plt

# 读取VASP计算产生的能带和态密度文件
vasprun = Vasprun("vasprun.xml")
band_structure = vasprun.get_band_structure("KPOINTS", line_mode=True, efermi=vasprun.efermi)
poscar = Poscar.from_file("POSCAR")

# 绘制能带图
plotter = BSPlotter(band_structure)
plotter.get_plot().show()
plt.savefig("band.png")
# 绘制态密度图
tdos = vasprun.tdos
dos = vasprun.complete_dos
proj_dos = band_structure.get_projection_on_elements()
plotter = DosPlotter(sigma=0.2)
plotter.add_dos("Total DOS", tdos)
plotter.add_dos_dict(proj_dos)
plotter.get_plot(xlim=[-10, 10]).show()
plt.savefig("dos.png")
