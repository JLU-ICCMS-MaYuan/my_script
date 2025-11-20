"""
自定义异常类模块

定义QE计算中可能出现的各种异常类型，便于精确的错误处理和诊断。
所有异常都继承自QEException基类。

作者：Claude
创建时间：2025-11-19
"""

from typing import Optional, List


# ============================================================================
# 基础异常类
# ============================================================================

class QEException(Exception):
    """QE模块的基础异常类

    所有自定义异常都应该继承自这个类，方便统一捕获和处理。

    Attributes
    ----------
    message : str
        错误信息
    details : dict, optional
        额外的错误详情（如文件路径、参数值等）
    """

    def __init__(self, message: str, details: Optional[dict] = None):
        self.message = message
        self.details = details or {}
        super().__init__(self.format_message())

    def format_message(self) -> str:
        """格式化错误信息"""
        msg = f"[QE Error] {self.message}"
        if self.details:
            msg += f"\n详细信息: {self.details}"
        return msg


# ============================================================================
# 文件相关异常
# ============================================================================

class FileNotFoundError(QEException):
    """文件未找到异常"""

    def __init__(self, file_path: str, file_type: Optional[str] = None):
        message = f"无法找到文件: {file_path}"
        if file_type:
            message += f" (文件类型: {file_type})"
        details = {'file_path': file_path, 'file_type': file_type}
        super().__init__(message, details)


class FileFormatError(QEException):
    """文件格式错误异常"""

    def __init__(self, file_path: str, expected_format: str, reason: Optional[str] = None):
        message = f"文件格式错误: {file_path}"
        message += f"\n期望格式: {expected_format}"
        if reason:
            message += f"\n错误原因: {reason}"
        details = {
            'file_path': file_path,
            'expected_format': expected_format,
            'reason': reason
        }
        super().__init__(message, details)


class FileParseError(QEException):
    """文件解析错误异常"""

    def __init__(self, file_path: str, line_number: Optional[int] = None, reason: Optional[str] = None):
        message = f"无法解析文件: {file_path}"
        if line_number:
            message += f"\n错误位置: 第 {line_number} 行"
        if reason:
            message += f"\n错误原因: {reason}"
        details = {
            'file_path': file_path,
            'line_number': line_number,
            'reason': reason
        }
        super().__init__(message, details)


# ============================================================================
# 配置相关异常
# ============================================================================

class ConfigurationError(QEException):
    """配置错误异常"""

    def __init__(self, parameter: str, value, reason: str):
        message = f"配置参数错误: {parameter} = {value}"
        message += f"\n错误原因: {reason}"
        details = {'parameter': parameter, 'value': value, 'reason': reason}
        super().__init__(message, details)


class MissingParameterError(QEException):
    """缺少必需参数异常"""

    def __init__(self, parameter: str, context: Optional[str] = None):
        message = f"缺少必需参数: {parameter}"
        if context:
            message += f"\n使用场景: {context}"
        details = {'parameter': parameter, 'context': context}
        super().__init__(message, details)


class InvalidParameterError(QEException):
    """无效参数异常"""

    def __init__(self, parameter: str, value, valid_range: Optional[str] = None, reason: Optional[str] = None):
        message = f"无效参数: {parameter} = {value}"
        if valid_range:
            message += f"\n有效范围: {valid_range}"
        if reason:
            message += f"\n错误原因: {reason}"
        details = {
            'parameter': parameter,
            'value': value,
            'valid_range': valid_range,
            'reason': reason
        }
        super().__init__(message, details)


# ============================================================================
# 计算相关异常
# ============================================================================

class CalculationError(QEException):
    """计算错误异常"""

    def __init__(self, calculation_type: str, reason: str, output_file: Optional[str] = None):
        message = f"{calculation_type} 计算失败"
        message += f"\n错误原因: {reason}"
        if output_file:
            message += f"\n输出文件: {output_file}"
        details = {
            'calculation_type': calculation_type,
            'reason': reason,
            'output_file': output_file
        }
        super().__init__(message, details)


class ConvergenceError(QEException):
    """收敛性错误异常"""

    def __init__(self, calculation_type: str, max_iterations: int, current_value: Optional[float] = None, threshold: Optional[float] = None):
        message = f"{calculation_type} 未收敛"
        message += f"\n最大迭代次数: {max_iterations}"
        if current_value and threshold:
            message += f"\n当前值: {current_value:.6e}, 阈值: {threshold:.6e}"
        details = {
            'calculation_type': calculation_type,
            'max_iterations': max_iterations,
            'current_value': current_value,
            'threshold': threshold
        }
        super().__init__(message, details)


class SCFNotConvergedError(ConvergenceError):
    """SCF自洽计算未收敛异常"""

    def __init__(self, max_iterations: int, final_accuracy: Optional[float] = None, threshold: Optional[float] = None):
        super().__init__(
            calculation_type="SCF自洽计算",
            max_iterations=max_iterations,
            current_value=final_accuracy,
            threshold=threshold
        )


class RelaxNotConvergedError(ConvergenceError):
    """结构优化未收敛异常"""

    def __init__(self, max_iterations: int, final_force: Optional[float] = None, force_threshold: Optional[float] = None):
        super().__init__(
            calculation_type="结构优化",
            max_iterations=max_iterations,
            current_value=final_force,
            threshold=force_threshold
        )


# ============================================================================
# 结构相关异常
# ============================================================================

class StructureError(QEException):
    """结构错误异常"""

    def __init__(self, reason: str, structure_file: Optional[str] = None):
        message = f"结构错误: {reason}"
        if structure_file:
            message += f"\n结构文件: {structure_file}"
        details = {'reason': reason, 'structure_file': structure_file}
        super().__init__(message, details)


class SymmetryError(QEException):
    """对称性错误异常"""

    def __init__(self, reason: str, expected_symmetry: Optional[str] = None, found_symmetry: Optional[str] = None):
        message = f"对称性错误: {reason}"
        if expected_symmetry:
            message += f"\n期望对称性: {expected_symmetry}"
        if found_symmetry:
            message += f"\n实际对称性: {found_symmetry}"
        details = {
            'reason': reason,
            'expected_symmetry': expected_symmetry,
            'found_symmetry': found_symmetry
        }
        super().__init__(message, details)


# ============================================================================
# 赝势相关异常
# ============================================================================

class PseudopotentialError(QEException):
    """赝势错误异常"""

    def __init__(self, element: str, reason: str, pseudopotential_dir: Optional[str] = None):
        message = f"元素 {element} 的赝势错误: {reason}"
        if pseudopotential_dir:
            message += f"\n赝势目录: {pseudopotential_dir}"
        details = {
            'element': element,
            'reason': reason,
            'pseudopotential_dir': pseudopotential_dir
        }
        super().__init__(message, details)


class PseudopotentialNotFoundError(PseudopotentialError):
    """赝势文件未找到异常"""

    def __init__(self, element: str, pseudopotential_dir: Optional[str] = None):
        super().__init__(
            element=element,
            reason="未找到合适的赝势文件",
            pseudopotential_dir=pseudopotential_dir
        )


# ============================================================================
# 声子计算相关异常
# ============================================================================

class PhononCalculationError(QEException):
    """声子计算错误异常"""

    def __init__(self, reason: str, q_point: Optional[List[float]] = None, mode: Optional[int] = None):
        message = f"声子计算错误: {reason}"
        if q_point:
            message += f"\nQ点: {q_point}"
        if mode:
            message += f"\n模式: {mode}"
        details = {'reason': reason, 'q_point': q_point, 'mode': mode}
        super().__init__(message, details)


class ImaginaryFrequencyError(QEException):
    """虚频错误异常"""

    def __init__(self, q_point: List[float], frequency: float, mode: int, threshold: float = -10.0):
        message = f"检测到虚频: {frequency:.2f} cm^-1"
        message += f"\nQ点: {q_point}, 模式: {mode}"
        message += f"\n建议: 检查结构优化是否充分，或增大degauss参数"
        details = {
            'q_point': q_point,
            'frequency': frequency,
            'mode': mode,
            'threshold': threshold
        }
        super().__init__(message, details)


# ============================================================================
# 工作流相关异常
# ============================================================================

class WorkflowError(QEException):
    """工作流错误异常"""

    def __init__(self, workflow_name: str, step: str, reason: str):
        message = f"工作流 '{workflow_name}' 在步骤 '{step}' 失败"
        message += f"\n错误原因: {reason}"
        details = {
            'workflow_name': workflow_name,
            'step': step,
            'reason': reason
        }
        super().__init__(message, details)


class DependencyError(QEException):
    """依赖错误异常"""

    def __init__(self, current_step: str, required_step: str, missing_file: Optional[str] = None):
        message = f"步骤 '{current_step}' 的依赖未满足"
        message += f"\n需要先完成步骤: {required_step}"
        if missing_file:
            message += f"\n缺少文件: {missing_file}"
        details = {
            'current_step': current_step,
            'required_step': required_step,
            'missing_file': missing_file
        }
        super().__init__(message, details)


# ============================================================================
# 任务调度相关异常
# ============================================================================

class SchedulerError(QEException):
    """任务调度错误异常"""

    def __init__(self, reason: str, task_id: Optional[str] = None):
        message = f"任务调度错误: {reason}"
        if task_id:
            message += f"\n任务ID: {task_id}"
        details = {'reason': reason, 'task_id': task_id}
        super().__init__(message, details)


class JobSubmissionError(QEException):
    """任务提交错误异常"""

    def __init__(self, job_system: str, reason: str, script_path: Optional[str] = None):
        message = f"无法提交任务到 {job_system}"
        message += f"\n错误原因: {reason}"
        if script_path:
            message += f"\n提交脚本: {script_path}"
        details = {
            'job_system': job_system,
            'reason': reason,
            'script_path': script_path
        }
        super().__init__(message, details)


# ============================================================================
# 绘图相关异常
# ============================================================================

class PlottingError(QEException):
    """绘图错误异常"""

    def __init__(self, plot_type: str, reason: str, data_file: Optional[str] = None):
        message = f"{plot_type} 绘图失败"
        message += f"\n错误原因: {reason}"
        if data_file:
            message += f"\n数据文件: {data_file}"
        details = {
            'plot_type': plot_type,
            'reason': reason,
            'data_file': data_file
        }
        super().__init__(message, details)


class DataNotFoundError(QEException):
    """数据未找到异常"""

    def __init__(self, data_type: str, expected_file: Optional[str] = None):
        message = f"未找到 {data_type} 数据"
        if expected_file:
            message += f"\n期望文件: {expected_file}"
        details = {'data_type': data_type, 'expected_file': expected_file}
        super().__init__(message, details)


# ============================================================================
# 辅助函数
# ============================================================================

def format_error_suggestion(error: Exception) -> str:
    """
    根据异常类型给出错误建议

    Parameters
    ----------
    error : Exception
        捕获的异常

    Returns
    -------
    str
        格式化的错误信息和建议
    """
    suggestions = {
        SCFNotConvergedError: [
            "1. 增加最大迭代次数 (electron_maxstep)",
            "2. 调整混合参数 (mixing_beta)",
            "3. 使用更密的k点网格",
            "4. 检查结构是否合理"
        ],
        RelaxNotConvergedError: [
            "1. 增加最大优化步数 (nstep)",
            "2. 调整力收敛阈值 (forc_conv_thr)",
            "3. 检查初始结构是否过于偏离平衡态",
            "4. 考虑使用更合适的优化算法"
        ],
        ImaginaryFrequencyError: [
            "1. 确保结构已充分优化（力小于0.001 Ry/Bohr）",
            "2. 增大展宽参数 degauss (如从0.02增加到0.05)",
            "3. 使用更密的k点网格进行SCF计算",
            "4. 如果是Gamma点虚频，可能结构不稳定，需重新优化"
        ],
        PseudopotentialNotFoundError: [
            "1. 检查赝势目录路径是否正确",
            "2. 确认该元素的赝势文件已下载",
            "3. 检查元素符号大小写是否正确",
            "4. 使用 get_USPP() 方法查看可用赝势"
        ]
    }

    error_type = type(error)
    suggestion_list = suggestions.get(error_type, [])

    msg = f"\n{'='*60}\n"
    msg += f"错误类型: {error_type.__name__}\n"
    msg += f"错误信息: {str(error)}\n"

    if suggestion_list:
        msg += f"\n建议的解决方案:\n"
        for suggestion in suggestion_list:
            msg += f"  {suggestion}\n"

    msg += f"{'='*60}\n"

    return msg


if __name__ == "__main__":
    # 测试示例
    print("QE异常类测试\n")

    # 测试1: 文件未找到
    try:
        raise FileNotFoundError("/path/to/structure.vasp", "POSCAR")
    except QEException as e:
        print(e)
        print()

    # 测试2: SCF未收敛
    try:
        raise SCFNotConvergedError(max_iterations=100, final_accuracy=1e-5, threshold=1e-6)
    except QEException as e:
        print(format_error_suggestion(e))

    # 测试3: 虚频错误
    try:
        raise ImaginaryFrequencyError(q_point=[0.0, 0.0, 0.0], frequency=-15.5, mode=3)
    except QEException as e:
        print(format_error_suggestion(e))

    # 测试4: 配置错误
    try:
        raise InvalidParameterError(
            parameter="ecutwfc",
            value=-10,
            valid_range="> 0",
            reason="截断能必须为正数"
        )
    except QEException as e:
        print(e)
