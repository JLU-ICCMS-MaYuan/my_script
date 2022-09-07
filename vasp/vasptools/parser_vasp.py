#!/usr/bin/env python

import numpy as np

from vasp.vasptools import poscar
from vasp.vasptools import vasprunxml
from vasp.vasptools import outcar
from structuregenerator.psolib.utils.uniq_atom_names import uniq_atom_names


class VASPPoscarFormat:
    def __init__(self):
        pass

    def from_poscar(self, file_name, **kwargs):
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
        data = poscar.to_system_data(lines)
        data = uniq_atom_names(data)
        return data

    def to_poscar(self, data, file_name, frame_idx=0, **kwargs):
        """
        Dump the system in vasp POSCAR format
        Parameters
        ----------
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        """
        w_str = VASPStringFormat().to_poscar(data, frame_idx=frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)


class VASPStringFormat:
    def __init__(self):
        pass

    def to_poscar(self, data, frame_idx=0, **kwargs):
        """
        Dump the system in vasp POSCAR format string
        Parameters
        ----------
        frame_idx : int
            The index of the frame to dump
        """
        assert frame_idx < len(data['coords'])
        return poscar.from_system_data(data, frame_idx)


# rotate the system to lammps convention
class VASPOutcarFormat:
    def __init__(self):
        pass

    def from_outcar(self, file_name, begin=0, step=1, **kwargs):
        data = {}
        ml = kwargs.get("ml", False)
        (
            data['atom_names'],
            data['atom_numbs'],
            data['atom_types'],
            data['cells'],
            data['coords'],
            data['energies'],
            data['forces'],
            tmp_virial,
        ) = outcar.get_frames(file_name, begin=begin, step=step, ml=ml)
        if tmp_virial is not None:
            data['virials'] = tmp_virial
        # scale virial to the unit of eV
        if 'virials' in data:
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data['cells'].shape[0]):
                vol = np.linalg.det(np.reshape(data['cells'][ii], [3, 3]))
                data['virials'][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        return data


# rotate the system to lammps convention
class VASPXMLFormat:
    def __init__(self):
        pass

    def from_vasprun_xml(self, file_name, begin=0, step=1, **kwargs):
        data = {}
        (
            data['atom_names'],
            data['atom_types'],
            data['cells'],
            data['coords'],
            data['energies'],
            data['forces'],
            data['virials'],
        ) = vasprunxml.analyze(file_name, type_idx_zero=True, begin=begin, step=step)
        data['atom_numbs'] = []
        for ii in range(len(data['atom_names'])):
            data['atom_numbs'].append(sum(data['atom_types'] == ii))
        # the vasp xml assumes the direct coordinates
        # apply the transform to the cartesan coordinates
        for ii in range(data['cells'].shape[0]):
            data['coords'][ii] = np.matmul(data['coords'][ii], data['cells'][ii])
        # scale virial to the unit of eV
        v_pref = 1 * 1e3 / 1.602176621e6
        for ii in range(data['cells'].shape[0]):
            vol = np.linalg.det(np.reshape(data['cells'][ii], [3, 3]))
            data['virials'][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        return data
