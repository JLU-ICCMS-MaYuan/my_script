"""
单元测试：常数和单位转换模块

作者：Claude
创建时间：2025-11-19
"""

import pytest
from core.constants import EnergyConversion, LengthConversion


class TestEnergyConversion:
    """测试能量单位转换"""

    def test_ry_to_ev(self):
        """测试Rydberg到eV的转换"""
        result = EnergyConversion.convert(1.0, 'Ry', 'eV')
        assert abs(result - 13.605693) < 0.0001

    def test_ev_to_ry(self):
        """测试eV到Rydberg的转换"""
        result = EnergyConversion.convert(13.605693, 'eV', 'Ry')
        assert abs(result - 1.0) < 0.0001

    def test_thz_to_cm_inv(self):
        """测试THz到cm^-1的转换"""
        result = EnergyConversion.convert(1.0, 'THz', 'cm-1')
        assert abs(result - 33.356410) < 0.001


class TestLengthConversion:
    """测试长度单位转换"""

    def test_bohr_to_angstrom(self):
        """测试Bohr到Angstrom的转换"""
        result = LengthConversion.convert(1.0, 'Bohr', 'Angstrom')
        assert abs(result - 0.529177) < 0.00001

    def test_angstrom_to_bohr(self):
        """测试Angstrom到Bohr的转换"""
        result = LengthConversion.convert(0.529177, 'Angstrom', 'Bohr')
        assert abs(result - 1.0) < 0.00001
