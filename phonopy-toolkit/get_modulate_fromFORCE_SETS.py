# """Tests for Modulation."""
from phonopy import Phonopy
import phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.cells import get_supercell 


# symbols_LaBH8 = ["La"] * 1 + ["B"] * 1 + ["H"] * 8
# cell_LaBH8 = PhonopyAtoms(
    # cell=lattice_LaBH8, scaled_positions=positions_LaBH8, symbols=symbols_LaBH8
# )

# phonon = phonopy.load(
#     supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
#     primitive_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#     unitcell_filename="POSCAR-init",
#     force_sets_filename="FORCE_SETS",
#     # born_filename="BORN",
# )
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS


unitcell = read_vasp("POSCAR-init")
phonon = Phonopy(
    unitcell=unitcell,
    supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    primitive_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
)
force_sets = parse_FORCE_SETS(filename="FORCE_SETS")
phonon.set_displacement_dataset(force_sets)
phonon.produce_force_constants()

# smat = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
# supercell_LaBH8 = get_supercell(cell_LaBH8, smat, is_old_style=True)

# supercell_LaBH8 = Phonopy(
#     # cell=supercell_LaBH8.cell, scaled_positions=supercell_LaBH8.positions, symbols=supercell_LaBH8.symbols
#     supercell_LaBH8
# )
# print(supercell_LaBH8)

def test_modulation(ph_input: Phonopy): #, helper_methods):
    """Test to calculate modulation by LaBH8."""
    ph = ph_input
    ph.run_modulations(
        dimension=[2, 2, 2], 
        phonon_modes=[[[0.5, 0.5, 0], 2, 2, 0], [[0.5, 0.5, 0], 2, 4, 0], [[0.5, 0.5, 0], 2, 6, 0]]
        # phonon_modes = [q-point,  band index(int),  amplitude(float),  phase(float)]
        )     
    cells = ph.get_modulated_supercells()
    print(cells[0].get_cell())
    # _show(cells)
    # helper_methods.compare_cells_with_order(cells[0], cell_LaBH8)


def _show(cells):
    for cell in cells:
        for p in cell.scaled_positions:
            print("[%.8f, %.8f, %.8f]," % tuple(p))


test_modulation(phonon)