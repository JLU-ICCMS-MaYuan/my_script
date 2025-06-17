import logging
import re

logger = logging.getLogger(__name__)

class qe_check:
    def __init__(self, output_file, check_file_type):
        """
        初始化 QECheck 类。
        
        参数:
            output_file (str): QE 输出文件的路径。
        """
        self.output_file = output_file
        self.content = self._read_output_file()
        
        self.check_file_type = check_file_type
        self.check_results = False
        if self.check_file_type == "scf":
            self.check_results = self.check_scf()
        elif self.check_file_type == "relax":
            self.check_results = self.check_relax()
            
    def _read_output_file(self):
        """
        读取 QE 输出文件的内容。
        
        返回:
            str: 输出文件的内容。
        """
        try:
            with open(self.output_file, 'r') as f:
                return f.read()
        except FileNotFoundError:
            logger.warning(f"file {self.output_file} 未找到。")
            return None
        
    def check_scf(self):
        """
        检查 SCF 计算是否完成。
        
        返回:
            bool: 如果 SCF 计算完成，返回 True; 否则返回 False。
        """
        if self.content is None:
            return False
        
        # 检查 SCF 是否完成的标志
        scf_completed = re.search(r"convergence has been achieved", self.content, re.IGNORECASE)
        if scf_completed:
            logger.debug("SCF 计算完成。")
            return True
        else:
            logger.debug("SCF 计算未完成。")
            return False

    def check_relax(self):
        """
        检查 RELAX 计算是否完成。
        
        返回:
            bool: 如果 RELAX 计算完成，返回 True；否则返回 False。
        """
        if self.content is None:
            return False

        # 检查 RELAX 是否完成的标志
        relax_completed = re.search(r"End final coordinates", self.content, re.IGNORECASE)
        if relax_completed:
            logger.debug("RELAX 计算完成。")
            return True
        else:
            logger.debug("RELAX 计算未完成。")
            return False

    def check_convergence(self):
        """
        检查计算是否收敛。
        
        返回:
            bool: 如果计算收敛，返回 True；否则返回 False。
        """
        if self.content is None:
            return False

        # 检查收敛标志
        converged = re.search(r"convergence has been achieved", self.content, re.IGNORECASE)
        if converged:
            logger.debug("计算收敛。")
            return True
        else:
            logger.debug("计算未收敛。")
            return False


# 示例用法
if __name__ == "__main__":
    output_file = "qe_output.out"  # 替换为你的 QE 输出文件路径
    qe_checker = qe_check(output_file)

    logger.debug("检查 SCF 计算：")
    qe_checker.check_scf()

    logger.debug("\n检查 RELAX 计算：")
    qe_checker.check_relax()

    logger.debug("\n检查计算是否收敛：")
    qe_checker.check_convergence()
