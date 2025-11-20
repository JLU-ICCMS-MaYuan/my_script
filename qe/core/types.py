"""
类型定义模块

定义QE计算中常用的类型别名和枚举类型，提高代码可读性和类型安全性。

作者：Claude
创建时间：2025-11-19
"""

from typing import List, Tuple, Dict, Union, Optional, Literal
from enum import Enum, auto
from pathlib import Path
from dataclasses import dataclass
import numpy as np


# ============================================================================
# 基础类型别名
# ============================================================================

# 路径类型（支持str和Path）
PathLike = Union[str, Path]

# 3D向量类型（晶胞向量、原子坐标等）
Vector3D = Union[List[float], Tuple[float, float, float], np.ndarray]

# 3x3矩阵类型（晶胞矩阵、旋转矩阵等）
Matrix3x3 = Union[List[List[float]], np.ndarray]

# K点网格类型
KPointGrid = Union[Tuple[int, int, int], List[int]]

# Q点网格类型
QPointGrid = Union[Tuple[int, int, int], List[int]]

# 元素组成类型 {元素名: 原子数}
Composition = Dict[str, int]


# ============================================================================
# 枚举类型
# ============================================================================

class CalculationType(str, Enum):
    """计算类型枚举"""
    RELAX = "relax"              # 固定晶胞结构优化
    RELAX_VC = "relax-vc"        # 可变晶胞结构优化
    SCF = "scf"                  # 自洽计算
    NSCF = "nscf"                # 非自洽计算
    BANDS = "bands"              # 能带计算
    PHONON = "ph"                # 声子计算
    DOS = "dos"                  # 态密度计算


class CoordinateType(str, Enum):
    """坐标类型枚举"""
    ALAT = "alat"                # 以晶格常数为单位
    BOHR = "bohr"                # 玻尔半径
    ANGSTROM = "angstrom"        # 埃
    CRYSTAL = "crystal"          # 晶体坐标（分数坐标）
    CARTESIAN = "cartesian"      # 笛卡尔坐标


class SmearingMethod(Enum):
    """展宽方法枚举"""
    GAUSSIAN = "gaussian"        # 高斯展宽
    METHFESSEL_PAXTON = "m-p"    # Methfessel-Paxton方法
    MARZARI_VANDERBILT = "m-v"   # Marzari-Vanderbilt冷涂抹
    FERMI_DIRAC = "f-d"          # Fermi-Dirac统计
    FIXED = "fixed"              # 固定占据
    TETRAHEDRA = "tetrahedra"    # 四面体方法


class JobSystem(str, Enum):
    """作业调度系统枚举"""
    SLURM = "slurm"
    PBS = "pbs"
    LSF = "lsf"
    BASH = "bash"                # 本地运行


class TaskStatus(str, Enum):
    """任务状态枚举"""
    PENDING = "pending"          # 等待中
    RUNNING = "running"          # 运行中
    SUCCESS = "success"          # 成功完成
    FAILED = "failed"            # 失败
    CANCELLED = "cancelled"      # 已取消
    TIMEOUT = "timeout"          # 超时


class FileFormat(str, Enum):
    """文件格式枚举"""
    VASP = "vasp"                # VASP格式（POSCAR/CONTCAR）
    CIF = "cif"                  # 晶体信息文件
    XYZ = "xyz"                  # XYZ坐标
    PDB = "pdb"                  # 蛋白质数据库格式
    QE_INPUT = "qe-in"           # QE输入文件
    QE_OUTPUT = "qe-out"         # QE输出文件


class PhononMode(str, Enum):
    """声子计算模式枚举"""
    NOSPLIT = "nosplit"          # 不分q点计算
    SPLIT_DYN0 = "split_dyn0"    # 按dyn0拆分
    SPLIT_ASSIGN_Q = "split_assignQ"  # 按Q点分配拆分
    MERGE = "merge"              # 合并动力学矩阵
    Q2R = "q2r"                  # 实空间力常数
    MATDYN = "matdyn"            # 声子谱计算
    PHONODOS = "phonodos"        # 声子态密度


class PipelineType(str, Enum):
    """Pipeline流程类型枚举"""
    RELAX_PHONO_SC = "relax-phono-sc"        # 优化→声子→超导
    RELAX_ELECTRON = "relax-electron"        # 优化→电子性质
    FULL_ANALYSIS = "full-analysis"          # 完整分析
    CONVERGENCE_TEST = "convergence"         # 收敛性测试


class PlotType(str, Enum):
    """绘图类型枚举"""
    PHONON_BAND = "phonon-band"              # 声子能带
    PHONON_DOS = "phonon-dos"                # 声子态密度
    ELECTRON_BAND = "electron-band"          # 电子能带
    ELECTRON_DOS = "electron-dos"            # 电子态密度
    PROJECTED_BAND = "projected-band"        # 投影能带
    PROJECTED_DOS = "projected-dos"          # 投影态密度
    ALPHA2F = "alpha2f"                      # Alpha2F函数
    TC_CURVE = "tc-curve"                    # Tc温度曲线


# ============================================================================
# 数据类（使用dataclass简化结构）
# ============================================================================

@dataclass
class KPointsConfig:
    """K点配置

    Attributes
    ----------
    grid : KPointGrid
        K点网格，例如 (6, 6, 6)
    shift : Tuple[float, float, float], optional
        K点网格偏移，默认 (0, 0, 0)
    """
    grid: KPointGrid
    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self):
        """验证K点配置"""
        if len(self.grid) != 3:
            raise ValueError(f"K点网格必须是3个整数: {self.grid}")
        if any(k <= 0 for k in self.grid):
            raise ValueError(f"K点网格必须都为正数: {self.grid}")


@dataclass
class Cell:
    """晶胞信息

    Attributes
    ----------
    lattice_vectors : Matrix3x3
        晶格向量 [a1, a2, a3]，每行一个向量
    positions : List[Vector3D]
        原子位置列表
    species : List[str]
        原子种类列表（与positions对应）
    coordinate_type : CoordinateType
        坐标类型
    """
    lattice_vectors: Matrix3x3
    positions: List[Vector3D]
    species: List[str]
    coordinate_type: CoordinateType = CoordinateType.CRYSTAL

    def __post_init__(self):
        """验证晶胞信息"""
        if len(self.positions) != len(self.species):
            raise ValueError(
                f"原子位置数量({len(self.positions)})与原子种类数量({len(self.species)})不匹配"
            )


@dataclass
class QPointInfo:
    """Q点信息

    Attributes
    ----------
    coordinates : Vector3D
        Q点坐标
    weight : float
        权重
    index : int
        索引号
    """
    coordinates: Vector3D
    weight: float
    index: int


@dataclass
class PhononModeInfo:
    """声子模式信息

    Attributes
    ----------
    mode_number : int
        模式编号
    frequency : float
        频率 [cm^-1]
    gamma : float, optional
        线宽 [THz]
    irrep : str, optional
        不可约表示
    """
    mode_number: int
    frequency: float
    gamma: Optional[float] = None
    irrep: Optional[str] = None


@dataclass
class ConvergenceTestResult:
    """收敛性测试结果

    Attributes
    ----------
    parameter_name : str
        测试的参数名称（如 'ecutwfc', 'kpoints'）
    values : List[float]
        测试的参数值列表
    energies : List[float]
        对应的能量值列表
    converged_value : Optional[float]
        收敛值（如果找到）
    threshold : float
        收敛阈值
    """
    parameter_name: str
    values: List[float]
    energies: List[float]
    converged_value: Optional[float] = None
    threshold: float = 0.001  # 默认1 meV


@dataclass
class SuperconductivityData:
    """超导性质数据

    Attributes
    ----------
    lambda_value : float
        电声耦合常数 λ
    omega_log : float
        对数平均声子频率 ω_log [K]
    tc_mcmillan : Optional[float]
        McMillan公式计算的Tc [K]
    tc_allen_dynes : Optional[float]
        Allen-Dynes公式计算的Tc [K]
    tc_eliashberg : Optional[float]
        Eliashberg方程计算的Tc [K]
    mu : float
        屏蔽常数 μ*
    """
    lambda_value: float
    omega_log: float
    tc_mcmillan: Optional[float] = None
    tc_allen_dynes: Optional[float] = None
    tc_eliashberg: Optional[float] = None
    mu: float = 0.10


# ============================================================================
# 类型守卫函数
# ============================================================================

def is_valid_kpoint_grid(grid) -> bool:
    """检查是否是有效的K点网格"""
    try:
        if len(grid) != 3:
            return False
        return all(isinstance(k, int) and k > 0 for k in grid)
    except TypeError:
        return False


def is_valid_vector3d(vec) -> bool:
    """检查是否是有效的3D向量"""
    try:
        if len(vec) != 3:
            return False
        return all(isinstance(x, (int, float)) for x in vec)
    except TypeError:
        return False


def is_valid_matrix3x3(mat) -> bool:
    """检查是否是有效的3x3矩阵"""
    try:
        if len(mat) != 3:
            return False
        return all(len(row) == 3 and all(isinstance(x, (int, float)) for x in row) for row in mat)
    except TypeError:
        return False


# ============================================================================
# 类型转换函数
# ============================================================================

def to_path(path_like: PathLike) -> Path:
    """将路径字符串或Path对象转换为Path对象"""
    return Path(path_like)


def to_vector3d(vec: Union[List, Tuple, np.ndarray]) -> Tuple[float, float, float]:
    """将向量转换为标准的3元组格式"""
    if not is_valid_vector3d(vec):
        raise ValueError(f"无效的3D向量: {vec}")
    return tuple(float(x) for x in vec)


def to_kpoint_grid(grid: Union[str, List, Tuple]) -> Tuple[int, int, int]:
    """
    将K点网格转换为标准的3元组格式

    支持格式:
    - "6 6 6" (字符串)
    - [6, 6, 6] (列表)
    - (6, 6, 6) (元组)
    """
    if isinstance(grid, str):
        grid = [int(x) for x in grid.strip().split()]

    if not is_valid_kpoint_grid(grid):
        raise ValueError(f"无效的K点网格: {grid}")

    return tuple(int(k) for k in grid)


# ============================================================================
# 单位类型（用于显式标注物理量的单位）
# ============================================================================

@dataclass
class Energy:
    """能量物理量（带单位）"""
    value: float
    unit: Literal["eV", "Ry", "Ha", "meV"] = "eV"

    def to(self, target_unit: str) -> "Energy":
        """转换单位"""
        try:
            from .constants import EnergyConversion
        except ImportError:
            from constants import EnergyConversion
        converted_value = EnergyConversion.convert(self.value, self.unit, target_unit)
        return Energy(converted_value, target_unit)


@dataclass
class Length:
    """长度物理量（带单位）"""
    value: float
    unit: Literal["Angstrom", "Bohr", "nm"] = "Angstrom"

    def to(self, target_unit: str) -> "Length":
        """转换单位"""
        try:
            from .constants import LengthConversion
        except ImportError:
            from constants import LengthConversion
        converted_value = LengthConversion.convert(self.value, self.unit, target_unit)
        return Length(converted_value, target_unit)


@dataclass
class Temperature:
    """温度物理量（带单位）"""
    value: float
    unit: Literal["K", "C", "eV"] = "K"

    def to(self, target_unit: str) -> "Temperature":
        """转换单位"""
        try:
            from .constants import TemperatureConversion
        except ImportError:
            from constants import TemperatureConversion

        # K → C
        if self.unit == "K" and target_unit == "C":
            return Temperature(TemperatureConversion.kelvin_to_celsius(self.value), "C")
        # C → K
        elif self.unit == "C" and target_unit == "K":
            return Temperature(TemperatureConversion.celsius_to_kelvin(self.value), "K")
        # K → eV
        elif self.unit == "K" and target_unit == "eV":
            return Temperature(TemperatureConversion.kelvin_to_ev(self.value), "eV")
        # eV → K
        elif self.unit == "eV" and target_unit == "K":
            return Temperature(TemperatureConversion.ev_to_kelvin(self.value), "K")
        else:
            return self


if __name__ == "__main__":
    # 测试示例
    print("QE类型定义模块测试\n")

    # 测试1: K点配置
    print("【测试K点配置】")
    kpts = KPointsConfig(grid=(6, 6, 6), shift=(0.5, 0.5, 0.5))
    print(f"K点网格: {kpts.grid}")
    print(f"偏移: {kpts.shift}\n")

    # 测试2: 能量物理量
    print("【测试能量单位转换】")
    e1 = Energy(13.6, "eV")
    e2 = e1.to("Ry")
    print(f"{e1.value} {e1.unit} = {e2.value:.6f} {e2.unit}\n")

    # 测试3: 长度物理量
    print("【测试长度单位转换】")
    l1 = Length(1.0, "Bohr")
    l2 = l1.to("Angstrom")
    print(f"{l1.value} {l1.unit} = {l2.value:.6f} {l2.unit}\n")

    # 测试4: 温度物理量
    print("【测试温度单位转换】")
    t1 = Temperature(300, "K")
    t2 = t1.to("eV")
    print(f"{t1.value} {t1.unit} = {t2.value:.6e} {t2.unit}\n")

    # 测试5: 类型守卫
    print("【测试类型守卫】")
    print(f"is_valid_kpoint_grid([6, 6, 6]): {is_valid_kpoint_grid([6, 6, 6])}")
    print(f"is_valid_kpoint_grid([6, 6]): {is_valid_kpoint_grid([6, 6])}")
    print(f"is_valid_vector3d([1.0, 2.0, 3.0]): {is_valid_vector3d([1.0, 2.0, 3.0])}\n")

    # 测试6: 枚举类型
    print("【测试枚举类型】")
    print(f"计算类型: {CalculationType.RELAX.value}")
    print(f"作业系统: {JobSystem.SLURM.value}")
    print(f"任务状态: {TaskStatus.SUCCESS.value}")
