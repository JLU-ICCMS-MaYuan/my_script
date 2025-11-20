"""
物理常数和单位转换模块

包含QE计算中常用的物理常数和单位转换因子。
所有常数均采用国际标准值（CODATA 2018）。

作者：Claude
创建时间：2025-11-19
"""

from typing import Dict
import math

# ============================================================================
# 基本物理常数
# ============================================================================

# 普朗克常数 [J·s]
PLANCK_CONSTANT = 6.62607015e-34

# 约化普朗克常数 ℏ [J·s]
HBAR = PLANCK_CONSTANT / (2 * math.pi)

# 玻尔兹曼常数 [J/K]
BOLTZMANN_CONSTANT = 1.380649e-23

# 电子电荷 [C]
ELECTRON_CHARGE = 1.602176634e-19

# 光速 [m/s]
SPEED_OF_LIGHT = 299792458

# 电子质量 [kg]
ELECTRON_MASS = 9.1093837015e-31

# 原子质量单位 [kg]
ATOMIC_MASS_UNIT = 1.66053906660e-27

# 阿伏伽德罗常数 [mol^-1]
AVOGADRO_CONSTANT = 6.02214076e23

# ============================================================================
# QE专用常数
# ============================================================================

# Hartree能量 [eV]
HARTREE_TO_EV = 27.211386245988

# Rydberg能量 [eV]
RYDBERG_TO_EV = 13.605693122994

# Bohr半径 [Angstrom]
BOHR_TO_ANGSTROM = 0.529177210903

# ============================================================================
# 能量单位转换因子
# ============================================================================

class EnergyConversion:
    """能量单位转换

    支持的单位：
    - eV (电子伏特)
    - Ry (Rydberg)
    - Ha (Hartree)
    - J (焦耳)
    - meV (毫电子伏特)
    - THz (太赫兹，用于频率)
    - cm-1 (波数)
    """

    # 以eV为基准的转换因子
    EV_TO_RY = 1.0 / RYDBERG_TO_EV  # eV → Ry
    EV_TO_HA = 1.0 / HARTREE_TO_EV  # eV → Ha
    EV_TO_J = ELECTRON_CHARGE       # eV → J
    EV_TO_MEV = 1000.0              # eV → meV

    RY_TO_EV = RYDBERG_TO_EV        # Ry → eV
    RY_TO_HA = 0.5                  # Ry → Ha
    RY_TO_MEV = RYDBERG_TO_EV * 1000.0  # Ry → meV

    HA_TO_EV = HARTREE_TO_EV        # Ha → eV
    HA_TO_RY = 2.0                  # Ha → Ry
    HA_TO_MEV = HARTREE_TO_EV * 1000.0  # Ha → meV

    # 频率单位转换
    # ω (Ry) = E (Ry) / ℏ (Ry·s)
    # 1 Ry = 3.289841960508e15 Hz = 3289.841960508 THz
    RY_TO_THZ = 3289.841960508      # Ry → THz
    THZ_TO_RY = 1.0 / RY_TO_THZ     # THz → Ry

    # 波数单位转换
    # 1 THz = 33.356409519815 cm^-1
    THZ_TO_CM_INV = 33.356409519815 # THz → cm^-1
    CM_INV_TO_THZ = 1.0 / THZ_TO_CM_INV  # cm^-1 → THz
    RY_TO_CM_INV = RY_TO_THZ * THZ_TO_CM_INV  # Ry → cm^-1

    @staticmethod
    def convert(value: float, from_unit: str, to_unit: str) -> float:
        """
        通用能量单位转换

        Parameters
        ----------
        value : float
            要转换的数值
        from_unit : str
            源单位 {'eV', 'Ry', 'Ha', 'meV', 'THz', 'cm-1'}
        to_unit : str
            目标单位 {'eV', 'Ry', 'Ha', 'meV', 'THz', 'cm-1'}

        Returns
        -------
        float
            转换后的数值

        Examples
        --------
        >>> EnergyConversion.convert(1.0, 'Ry', 'eV')
        13.605693122994
        >>> EnergyConversion.convert(100.0, 'THz', 'cm-1')
        3335.6409519815
        """
        # 单位标准化
        unit_map = {
            'ev': 'eV', 'ry': 'Ry', 'ha': 'Ha',
            'mev': 'meV', 'thz': 'THz', 'cm-1': 'cm-1', 'cm^-1': 'cm-1'
        }
        from_unit = unit_map.get(from_unit.lower(), from_unit)
        to_unit = unit_map.get(to_unit.lower(), to_unit)

        if from_unit == to_unit:
            return value

        # 转换路径：所有单位 → eV → 目标单位
        conversion_to_ev = {
            'eV': 1.0,
            'Ry': EnergyConversion.RY_TO_EV,
            'Ha': EnergyConversion.HA_TO_EV,
            'meV': 1.0 / 1000.0,
            'THz': EnergyConversion.RY_TO_EV * EnergyConversion.THZ_TO_RY,
            'cm-1': EnergyConversion.RY_TO_EV * EnergyConversion.THZ_TO_RY * EnergyConversion.CM_INV_TO_THZ
        }

        conversion_from_ev = {
            'eV': 1.0,
            'Ry': EnergyConversion.EV_TO_RY,
            'Ha': EnergyConversion.EV_TO_HA,
            'meV': 1000.0,
            'THz': EnergyConversion.THZ_TO_RY * EnergyConversion.EV_TO_RY,
            'cm-1': EnergyConversion.CM_INV_TO_THZ * EnergyConversion.THZ_TO_RY * EnergyConversion.EV_TO_RY
        }

        # 先转换到eV
        value_in_ev = value * conversion_to_ev.get(from_unit, 1.0)
        # 再转换到目标单位
        return value_in_ev * conversion_from_ev.get(to_unit, 1.0)


# ============================================================================
# 长度单位转换因子
# ============================================================================

class LengthConversion:
    """长度单位转换

    支持的单位：
    - Angstrom (埃)
    - Bohr (玻尔半径)
    - nm (纳米)
    - m (米)
    """

    BOHR_TO_ANGSTROM = BOHR_TO_ANGSTROM  # Bohr → Å
    ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM  # Å → Bohr
    ANGSTROM_TO_NM = 0.1            # Å → nm
    NM_TO_ANGSTROM = 10.0           # nm → Å
    ANGSTROM_TO_M = 1e-10           # Å → m

    @staticmethod
    def convert(value: float, from_unit: str, to_unit: str) -> float:
        """
        通用长度单位转换

        Parameters
        ----------
        value : float
            要转换的数值
        from_unit : str
            源单位 {'Angstrom', 'Bohr', 'nm', 'm'}
        to_unit : str
            目标单位 {'Angstrom', 'Bohr', 'nm', 'm'}

        Returns
        -------
        float
            转换后的数值
        """
        # 单位标准化
        unit_map = {
            'angstrom': 'Angstrom', 'å': 'Angstrom', 'a': 'Angstrom',
            'bohr': 'Bohr', 'a.u.': 'Bohr',
            'nm': 'nm', 'nanometer': 'nm',
            'm': 'm', 'meter': 'm'
        }
        from_unit = unit_map.get(from_unit.lower(), from_unit)
        to_unit = unit_map.get(to_unit.lower(), to_unit)

        if from_unit == to_unit:
            return value

        # 转换路径：所有单位 → Angstrom → 目标单位
        to_angstrom = {
            'Angstrom': 1.0,
            'Bohr': LengthConversion.BOHR_TO_ANGSTROM,
            'nm': LengthConversion.NM_TO_ANGSTROM,
            'm': 1e10
        }

        from_angstrom = {
            'Angstrom': 1.0,
            'Bohr': LengthConversion.ANGSTROM_TO_BOHR,
            'nm': LengthConversion.ANGSTROM_TO_NM,
            'm': LengthConversion.ANGSTROM_TO_M
        }

        value_in_angstrom = value * to_angstrom.get(from_unit, 1.0)
        return value_in_angstrom * from_angstrom.get(to_unit, 1.0)


# ============================================================================
# 温度单位转换因子
# ============================================================================

class TemperatureConversion:
    """温度单位转换"""

    @staticmethod
    def kelvin_to_celsius(k: float) -> float:
        """K → °C"""
        return k - 273.15

    @staticmethod
    def celsius_to_kelvin(c: float) -> float:
        """°C → K"""
        return c + 273.15

    @staticmethod
    def kelvin_to_ev(k: float) -> float:
        """
        K → eV (通过 kB·T)

        用于将温度转换为能量单位
        """
        return k * BOLTZMANN_CONSTANT / ELECTRON_CHARGE

    @staticmethod
    def ev_to_kelvin(ev: float) -> float:
        """eV → K (通过 E / kB)"""
        return ev * ELECTRON_CHARGE / BOLTZMANN_CONSTANT


# ============================================================================
# 常用常数字典（方便查询）
# ============================================================================

PHYSICAL_CONSTANTS: Dict[str, Dict[str, float]] = {
    'Planck constant': {'value': PLANCK_CONSTANT, 'unit': 'J·s'},
    'reduced Planck constant': {'value': HBAR, 'unit': 'J·s'},
    'Boltzmann constant': {'value': BOLTZMANN_CONSTANT, 'unit': 'J/K'},
    'electron charge': {'value': ELECTRON_CHARGE, 'unit': 'C'},
    'speed of light': {'value': SPEED_OF_LIGHT, 'unit': 'm/s'},
    'electron mass': {'value': ELECTRON_MASS, 'unit': 'kg'},
    'atomic mass unit': {'value': ATOMIC_MASS_UNIT, 'unit': 'kg'},
    'Avogadro constant': {'value': AVOGADRO_CONSTANT, 'unit': 'mol^-1'},
    'Hartree energy': {'value': HARTREE_TO_EV, 'unit': 'eV'},
    'Rydberg energy': {'value': RYDBERG_TO_EV, 'unit': 'eV'},
    'Bohr radius': {'value': BOHR_TO_ANGSTROM, 'unit': 'Angstrom'},
}


# ============================================================================
# 辅助函数
# ============================================================================

def print_conversion_table():
    """打印常用单位转换表（用于参考）"""
    print("=" * 60)
    print("QE常用单位转换表")
    print("=" * 60)

    print("\n【能量单位】")
    print(f"1 Ry = {EnergyConversion.RY_TO_EV:.6f} eV")
    print(f"1 Ha = {EnergyConversion.HA_TO_EV:.6f} eV")
    print(f"1 eV = {EnergyConversion.EV_TO_RY:.6f} Ry")
    print(f"1 eV = {EnergyConversion.EV_TO_HA:.6f} Ha")

    print("\n【频率单位】")
    print(f"1 Ry = {EnergyConversion.RY_TO_THZ:.6f} THz")
    print(f"1 THz = {EnergyConversion.THZ_TO_CM_INV:.6f} cm^-1")
    print(f"1 Ry = {EnergyConversion.RY_TO_CM_INV:.6f} cm^-1")

    print("\n【长度单位】")
    print(f"1 Bohr = {LengthConversion.BOHR_TO_ANGSTROM:.6f} Angstrom")
    print(f"1 Angstrom = {LengthConversion.ANGSTROM_TO_BOHR:.6f} Bohr")
    print(f"1 Angstrom = {LengthConversion.ANGSTROM_TO_NM:.6f} nm")

    print("\n【温度-能量转换】")
    print(f"1 K = {TemperatureConversion.kelvin_to_ev(1.0):.6e} eV")
    print(f"300 K = {TemperatureConversion.kelvin_to_ev(300.0):.6f} eV")

    print("=" * 60)


if __name__ == "__main__":
    # 测试示例
    print("QE物理常数和单位转换模块测试\n")

    # 能量转换测试
    print("【能量转换测试】")
    print(f"1 Ry = {EnergyConversion.convert(1.0, 'Ry', 'eV'):.6f} eV")
    print(f"100 eV = {EnergyConversion.convert(100.0, 'eV', 'Ry'):.6f} Ry")
    print(f"1000 cm^-1 = {EnergyConversion.convert(1000.0, 'cm-1', 'THz'):.6f} THz")

    # 长度转换测试
    print("\n【长度转换测试】")
    print(f"1 Bohr = {LengthConversion.convert(1.0, 'Bohr', 'Angstrom'):.6f} Å")
    print(f"5 Å = {LengthConversion.convert(5.0, 'Angstrom', 'Bohr'):.6f} Bohr")

    # 打印转换表
    print("\n")
    print_conversion_table()
