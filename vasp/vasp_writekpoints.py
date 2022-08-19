import os

from pymatgen.io.vasp import Kpoints

from vasp_inputpara import vasp_inputpara


class vasp_writekpoints:
    def __init__(
        self, 
        work_underpressure: str,
        **kwargs, 
        ):
        self.work_underpressure = work_underpressure
        for key, value in kwargs.items():
            setattr(self, key, value)        

        kpoints = Kpoints.automatic_density_by_lengths(
            self.sposcar_struct_type, 
            length_densities=self.kpoints,
            force_gamma=True)
        kpoints.write_file(os.path.join(self.work_underpressure, "KPOINTS"))
        print(kpoints)
    
    @classmethod
    def init_from_inputpara(cls, other_class: vasp_inputpara):
        self = cls(
            work_underpressure=other_class.work_underpressure,
            sposcar_struct_type=other_class.sposcar_struct_type,
            kpoints=other_class.kpoints,
        )
        return self