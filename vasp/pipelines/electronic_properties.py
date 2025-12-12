"""
VASP电子性质Pipeline

完整流程：结构优化 → 自洽 → DOS → 能带 → ELF → COHP → 绘图

作者：Claude
创建时间：2025-11-20
"""

import logging
import shutil
from pathlib import Path
from typing import Optional, List

from vasp.pipelines.base import BasePipeline
from vasp.analysis import plotters
from vasp.utils.job import load_job_config, write_job_script, submit_job

logger = logging.getLogger(__name__)


class ElectronicPropertiesPipeline(BasePipeline):
    """
    电子性质全流程Pipeline

    执行步骤：
    1. 结构优化（Relax）
    2. 自洽计算（SCF）
    3. DOS计算
    4. 能带计算
    5. ELF计算（可选）
    6. COHP计算（可选）
    7. 自动绘图

    所有计算完全自动化，中间不需要人工干预
    """

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        kspacing: float = 0.2,
        encut: Optional[float] = None,
        include_elf: bool = True,
        include_cohp: bool = True,
        plot_dos_type: str = "element",
        queue_system: Optional[str] = None,
        mpi_procs: Optional[int] = None,
        potcar_dir: Optional[Path] = None,
        potcar_type: str = "PBE",
        **kwargs
    ):
        """
        初始化电子性质Pipeline

        Parameters
        ----------
        structure_file : Path
            输入结构文件（POSCAR格式）
        work_dir : Path
            工作目录
        kspacing : float
            K点间距（Angstrom^-1）
        encut : float, optional
            平面波截断能（eV），不指定则使用POTCAR推荐值
        include_elf : bool
            是否包含ELF计算
        include_cohp : bool
            是否包含COHP计算
        plot_dos_type : str
            DOS投影类型：'element', 'spd', 'element_spd'
        queue_system : str, optional
            队列系统：'slurm', 'pbs', 'bash'
        mpi_procs : int, optional
            mpirun -np 参数，默认取 rc 中的设置或8
        potcar_dir : Path, optional
            POTCAR库目录，如果不提供则需要手动准备POTCAR
        potcar_type : str
            POTCAR类型：'PBE', 'LDA', 'PW91'等
        """
        super().__init__(structure_file, work_dir, **kwargs)

        self.job_cfg = load_job_config()
        self.kspacing = kspacing
        self.encut = encut
        self.include_elf = include_elf
        self.include_cohp = include_cohp
        self.plot_dos_type = plot_dos_type
        self.queue_system = queue_system or "bash"
        self.mpi_procs = mpi_procs
        default_potcar = self.job_cfg.potcar_dir if self.job_cfg else None
        self.potcar_dir = Path(potcar_dir) if potcar_dir else default_potcar
        self.potcar_type = potcar_type

        # 子目录
        self.relax_dir = self.work_dir / "01_relax"
        self.scf_dir = self.work_dir / "02_scf"
        self.dos_dir = self.work_dir / "03_dos"
        self.band_dir = self.work_dir / "04_band"
        self.elf_dir = self.work_dir / "05_elf"
        self.cohp_dir = self.work_dir / "06_cohp"
        self.plots_dir = self.work_dir / "plots"

    def get_steps(self) -> List[str]:
        """返回所有步骤"""
        steps = [
            "relax",
            "scf",
            "dos",
            "band",
        ]

        if self.include_elf:
            steps.append("elf")

        if self.include_cohp:
            steps.append("cohp")

        steps.append("plotting")

        return steps

    def execute_step(self, step_name: str) -> bool:
        """执行单个步骤"""
        try:
            if step_name == "relax":
                return self._run_relax()
            elif step_name == "scf":
                return self._run_scf()
            elif step_name == "dos":
                return self._run_dos()
            elif step_name == "band":
                return self._run_band()
            elif step_name == "elf":
                return self._run_elf()
            elif step_name == "cohp":
                return self._run_cohp()
            elif step_name == "plotting":
                return self._run_plotting()
            else:
                logger.error(f"未知步骤: {step_name}")
                return False

        except Exception as e:
            logger.error(f"步骤 '{step_name}' 执行异常: {e}", exc_info=True)
            return False

    def _run_relax(self) -> bool:
        """Step 1: 结构优化"""
        logger.info("执行结构优化...")

        self.relax_dir.mkdir(parents=True, exist_ok=True)

        # 复制POSCAR
        shutil.copy(self.structure_file, self.relax_dir / "POSCAR")

        # 创建INCAR（结构优化）
        self._write_relax_incar(self.relax_dir / "INCAR")

        # 创建KPOINTS
        self._write_kpoints(self.relax_dir / "KPOINTS", self.kspacing)

        # 准备POTCAR
        if self.potcar_dir:
            from vasp.pipelines.utils import prepare_potcar
            if not prepare_potcar(
                self.relax_dir / "POSCAR",
                self.potcar_dir,
                self.relax_dir / "POTCAR",
                self.potcar_type
            ):
                logger.error("POTCAR准备失败")
                return False
        else:
            logger.warning("未提供potcar_dir，请确保POTCAR文件已手动准备好")

        # 提交任务
        job_script = self._write_job_script(self.relax_dir, "relax")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.relax_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 relax 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.relax_dir, self.queue_system):
            return False

        # 检查收敛
        if not self._check_convergence(self.relax_dir):
            logger.error("结构优化未收敛")
            return False

        # 保存优化后的结构
        contcar = self.relax_dir / "CONTCAR"
        if contcar.exists():
            shutil.copy(contcar, self.work_dir / "POSCAR_relaxed")
            self.steps_data['relaxed_structure'] = str(self.work_dir / "POSCAR_relaxed")

        logger.info("结构优化完成")
        return True

    def _run_scf(self) -> bool:
        """Step 2: 自洽计算"""
        logger.info("执行自洽计算...")

        self.scf_dir.mkdir(parents=True, exist_ok=True)

        # 使用优化后的结构
        relaxed_poscar = Path(self.steps_data.get('relaxed_structure', self.structure_file))
        shutil.copy(relaxed_poscar, self.scf_dir / "POSCAR")

        # 创建INCAR（高精度SCF）
        self._write_scf_incar(self.scf_dir / "INCAR")

        # 创建KPOINTS（更密）
        self._write_kpoints(self.scf_dir / "KPOINTS", self.kspacing)

        # 复制POTCAR（从relax目录）
        potcar_source = self.relax_dir / "POTCAR"
        if potcar_source.exists():
            shutil.copy(potcar_source, self.scf_dir / "POTCAR")
        else:
            logger.warning(f"未找到POTCAR文件: {potcar_source}")

        # 提交任务
        job_script = self._write_job_script(self.scf_dir, "scf")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.scf_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 scf 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.scf_dir, self.queue_system):
            return False

        # 保存CHGCAR供后续使用
        chgcar = self.scf_dir / "CHGCAR"
        if chgcar.exists():
            self.steps_data['chgcar'] = str(chgcar)

        logger.info("自洽计算完成")
        return True

    def _run_dos(self) -> bool:
        """Step 3: DOS计算"""
        logger.info("执行DOS计算...")

        self.dos_dir.mkdir(parents=True, exist_ok=True)

        # 复制文件
        shutil.copy(self.scf_dir / "POSCAR", self.dos_dir / "POSCAR")
        shutil.copy(self.scf_dir / "CHGCAR", self.dos_dir / "CHGCAR")
        shutil.copy(self.scf_dir / "POTCAR", self.dos_dir / "POTCAR")

        # 创建INCAR（DOS）
        self._write_dos_incar(self.dos_dir / "INCAR")

        # 创建KPOINTS（DOS需要更密的网格）
        self._write_kpoints(self.dos_dir / "KPOINTS", self.kspacing / 2)

        # 提交任务
        job_script = self._write_job_script(self.dos_dir, "dos")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.dos_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 dos 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.dos_dir, self.queue_system):
            return False

        logger.info("DOS计算完成")
        return True

    def _run_band(self) -> bool:
        """Step 4: 能带计算"""
        logger.info("执行能带计算...")

        self.band_dir.mkdir(parents=True, exist_ok=True)

        # 复制文件
        shutil.copy(self.scf_dir / "POSCAR", self.band_dir / "POSCAR")
        shutil.copy(self.scf_dir / "CHGCAR", self.band_dir / "CHGCAR")
        shutil.copy(self.scf_dir / "POTCAR", self.band_dir / "POTCAR")

        # 创建INCAR（能带）
        self._write_band_incar(self.band_dir / "INCAR")

        # 创建KPOINTS（高对称路径）
        self._write_band_kpoints(self.band_dir / "KPOINTS")

        # 提交任务
        job_script = self._write_job_script(self.band_dir, "band")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.band_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 band 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.band_dir, self.queue_system):
            return False

        logger.info("能带计算完成")
        return True

    def _run_elf(self) -> bool:
        """Step 5: ELF计算"""
        logger.info("执行ELF计算...")

        self.elf_dir.mkdir(parents=True, exist_ok=True)

        # 复制文件
        shutil.copy(self.scf_dir / "POSCAR", self.elf_dir / "POSCAR")
        shutil.copy(self.scf_dir / "CHGCAR", self.elf_dir / "CHGCAR")
        shutil.copy(self.scf_dir / "POTCAR", self.elf_dir / "POTCAR")

        # 创建INCAR（ELF）
        self._write_elf_incar(self.elf_dir / "INCAR")

        # KPOINTS使用SCF的
        shutil.copy(self.scf_dir / "KPOINTS", self.elf_dir / "KPOINTS")

        # 提交任务
        job_script = self._write_job_script(self.elf_dir, "elf")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.elf_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 elf 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.elf_dir, self.queue_system):
            return False

        logger.info("ELF计算完成")
        return True

    def _run_cohp(self) -> bool:
        """Step 6: COHP计算"""
        logger.info("执行COHP计算...")

        self.cohp_dir.mkdir(parents=True, exist_ok=True)

        # 复制文件
        shutil.copy(self.scf_dir / "POSCAR", self.cohp_dir / "POSCAR")
        shutil.copy(self.scf_dir / "POTCAR", self.cohp_dir / "POTCAR")

        # 创建INCAR（COHP）
        self._write_cohp_incar(self.cohp_dir / "INCAR")

        # KPOINTS
        shutil.copy(self.scf_dir / "KPOINTS", self.cohp_dir / "KPOINTS")

        # 提交任务
        job_script = self._write_job_script(self.cohp_dir, "cohp")
        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.cohp_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 cohp 作业，退出等待。")
            return True

        # 等待完成
        if not self._wait_for_job(job_id, self.cohp_dir, self.queue_system):
            return False

        logger.info("COHP计算完成")
        return True

    def _run_plotting(self) -> bool:
        """Step 7: 自动绘图"""
        logger.info("开始绘图...")

        self.plots_dir.mkdir(parents=True, exist_ok=True)

        try:
            # 绘制能带
            plotters.plot_band_structure(
                self.band_dir,
                self.plots_dir / "band.png"
            )

            # 绘制DOS
            plotters.plot_dos(
                self.dos_dir,
                self.plots_dir / "dos.png",
                pdos_type=self.plot_dos_type
            )

            # 绘制ELF（如果有）
            if self.include_elf and self.elf_dir.exists():
                try:
                    plotters.plot_elf(
                        self.elf_dir,
                        self.plots_dir / "elf.png"
                    )
                except Exception as e:
                    logger.warning(f"ELF绘图失败: {e}")

            # 绘制COHP（如果有）
            if self.include_cohp and self.cohp_dir.exists():
                try:
                    plotters.plot_cohp(
                        self.cohp_dir,
                        self.plots_dir / "cohp.png"
                    )
                except Exception as e:
                    logger.warning(f"COHP绘图失败: {e}")

            logger.info(f"所有图表已保存到: {self.plots_dir}")
            return True

        except Exception as e:
            logger.error(f"绘图失败: {e}", exc_info=True)
            return False

    def _write_relax_incar(self, incar_file: Path):
        """写入结构优化INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# Relaxation INCAR\n")
            f.write("SYSTEM = Structure Relaxation\n\n")
            f.write("# Electronic\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("EDIFF = 1E-6\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n\n")
            f.write("# Ionic\n")
            f.write("IBRION = 2\n")
            f.write("NSW = 200\n")
            f.write("ISIF = 3\n")
            f.write("EDIFFG = -0.01\n\n")
            f.write("# Output\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .TRUE.\n")

    def _write_scf_incar(self, incar_file: Path):
        """写入SCF INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# SCF INCAR\n")
            f.write("SYSTEM = Self-Consistent Calculation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("EDIFF = 1E-6\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .TRUE.\n")

    def _write_dos_incar(self, incar_file: Path):
        """写入DOS INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# DOS INCAR\n")
            f.write("SYSTEM = DOS Calculation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("ICHARG = 11\n")
            f.write("ISMEAR = -5\n")
            f.write("LORBIT = 11\n")
            f.write("NEDOS = 2000\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_band_incar(self, incar_file: Path):
        """写入能带INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# Band Structure INCAR\n")
            f.write("SYSTEM = Band Structure\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("ICHARG = 11\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n")
            f.write("LORBIT = 11\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_elf_incar(self, incar_file: Path):
        """写入ELF INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# ELF INCAR\n")
            f.write("SYSTEM = ELF Calculation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("ICHARG = 11\n")
            f.write("LELF = .TRUE.\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_cohp_incar(self, incar_file: Path):
        """写入COHP INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# COHP INCAR\n")
            f.write("SYSTEM = COHP Calculation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("ISMEAR = -5\n")
            f.write("LORBIT = 11\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_kpoints(self, kpoints_file: Path, kspacing: float):
        """写入自动K点"""
        with open(kpoints_file, 'w') as f:
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{kspacing} {kspacing} {kspacing}\n")

    def _write_band_kpoints(self, kpoints_file: Path):
        """写入能带K点（高对称路径）"""
        # 这里简化处理，使用seekpath或pymatgen自动生成
        # 实际使用时需要根据晶体结构确定
        with open(kpoints_file, 'w') as f:
            f.write("k-points for band structure\n")
            f.write("10\n")
            f.write("Line-mode\n")
            f.write("reciprocal\n")
            f.write("0.0 0.0 0.0   ! Gamma\n")
            f.write("0.5 0.0 0.0   ! X\n\n")
            f.write("0.5 0.0 0.0   ! X\n")
            f.write("0.5 0.5 0.0   ! M\n\n")
            f.write("0.5 0.5 0.0   ! M\n")
            f.write("0.0 0.0 0.0   ! Gamma\n")

    def _write_job_script(self, work_dir: Path, job_name: str) -> str:
        """写入任务提交脚本（支持 bash/slurm/pbs/lsf）"""
        script_path = write_job_script(
            work_dir=work_dir,
            job_name=job_name,
            queue_system=self.queue_system,
            cfg=self.job_cfg,
            use_gamma=False,
            mpi_procs=self.mpi_procs,
        )
        return str(script_path)

    def _submit_job(self, work_dir: Path, job_script: str) -> str:
        """提交任务"""
        return submit_job(Path(job_script), self.queue_system)
