#!/work/home/mayuan/miniconda3/envs/cage/bin/python3
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, PDPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element
import plotly.express as px
import pandas as pd

convexhull_data = pd.read_csv('nnconvexhull.csv', header=0, sep=',')

ini_entries = []
for idx, row in convexhull_data.iterrows():
    comp = Composition(row['formula'])
    enth = row['enthalpy']
    entry_id = row['Number']
    _entry = PDEntry(comp, enth)
    _entry.entry_id = entry_id
    ini_entries.append(_entry)

ini_pd = PhaseDiagram(ini_entries)
plotter = PDPlotter(ini_pd, show_unstable=0.2, backend='plotly')
plotter.show()
plotter.write_image('pd.png', image_format='png')
