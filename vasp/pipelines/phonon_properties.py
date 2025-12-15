"""
VASP声子性质Pipeline

完整流程：结构优化 → 声子谱计算 → 声子DOS → 绘图

作者：Claude
创建时间：2025-11-20
"""

import os
import logging
import shutil
from pathlib import Path
from typing import Optional, List

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import ensure_poscar, prepare_potcar, check_vasp_completion
from vasp.analysis import plotters
from vasp.utils.job import load_job_config, write_job_script, submit_job

logger = logging.getLogger(__name__)


class PhononPropertiesPipeline(BasePipeline):
    """
    声子性质全流程Pipeline

    执行步骤：
    1. 结构优化（Relax）
    2. 声子谱计算（使用位移法disp或DFPT法）
    3. Phonopy后处理（声子能带）
    4. Phonopy后处理（声子DOS）
    5. 自动绘图

    所有计算完全自动化
    """

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        supercell: List[int] = None,
        method: str = "disp",
        kdensity: Optional[float] = None,
        kspacing: Optional[float] = 0.3,
        encut: Optional[float] = None,
        queue_system: Optional[str] = None,
        mpi_procs: Optional[str] = None,
        potcar_dir: Optional[Path] = None,
        potcar_type: str = "PBE",
        include_relax: bool = True,
        custom_steps: Optional[List[str]] = None,
        pressure: float = 0.0,
        **kwargs
    ):
        """
        初始化声子性质Pipeline

        Parameters
        ----------
        structure_file : Path
            输入结构文件（POSCAR格式）
        work_dir : Path
            工作目录
        supercell : List[int]
            超胞大小，如[2, 2, 2]
        method : str
            声子计算方法：'disp'（位移法）或'dfpt'（密度泛函微扰理论）
        kdensity : float, optional
            K点密度（用于phonopy）
        kspacing : float
            K点间距
        encut : float, optional
            平面波截断能
        queue_system : str, optional
            队列系统
        mpi_procs : int, optional
            mpirun -np 参数，默认取 rc 中设置或8
        potcar_dir : Path, optional
            POTCAR库目录，如果不提供则需要手动准备POTCAR
        potcar_type : str
            POTCAR类型：'PBE', 'LDA', 'PW91'等
        include_relax : bool
            是否在声子前自动进行结构优化
        custom_steps : List[str], optional
            自定义步骤序列，如 ["phonon_prepare", "phonon_calculate"]
        pressure : float
            施加的外压（GPa），写入 PSTRESS（kBar）
        """
        super().__init__(structure_file, work_dir, pressure=pressure, **kwargs)

        self.job_cfg = load_job_config()
        self.supercell = supercell or [2, 2, 2]
        self.method = method
        self.kdensity = kdensity or 8000
        self.kspacing = kspacing
        self.encut = encut
        self.queue_system = queue_system or "bash"
        self.mpi_procs = mpi_procs
        # 声子必须基于结构优化
        self.include_relax = include_relax
        self.custom_steps = self._normalize_steps(custom_steps)
        self.pressure = pressure
        default_potcar = self.job_cfg.potcar_dir if self.job_cfg else None
        self.potcar_dir = Path(potcar_dir) if potcar_dir else default_potcar
        self.potcar_type = potcar_type

        # 子目录
        self.relax_dir = self.work_dir / "01_relax"
        self.phonon_dir = self.work_dir / "02_phonon"
        self.plots_dir = self.work_dir / "plots"

    def _normalize_steps(self, custom_steps: Optional[List[str]]) -> Optional[List[str]]:
        if not custom_steps:
            return None
        allowed = [
            "relax",
            "phonon_prepare",
            "phonon_calculate",
            "phonon_band",
            "phonon_dos",
            "plotting",
        ]
        normalized: List[str] = []
        for step in custom_steps:
            name = step.strip().lower()
            if name in allowed and name not in normalized:
                normalized.append(name)
            else:
                if name not in allowed:
                    logger.warning(f"忽略未支持的步骤: {name}")
        if normalized:
            if self.include_relax and "relax" not in normalized:
                normalized.insert(0, "relax")
            return normalized
        return None

    def get_steps(self) -> List[str]:
        """返回所有步骤"""
        if self.custom_steps:
            return self.custom_steps

        steps = [
            "phonon_prepare",
            "phonon_calculate",
            "phonon_band",
            "phonon_dos",
            "plotting",
        ]
        if self.include_relax:
            steps.insert(0, "relax")
        return steps

    def execute_step(self, step_name: str) -> bool:
        """执行单个步骤"""
        try:
            if step_name == "relax":
                return self._run_relax()
            elif step_name == "phonon_prepare":
                return self._prepare_phonon()
            elif step_name == "phonon_calculate":
                return self._run_phonon_calculate()
            elif step_name == "phonon_band":
                return self._run_phonon_band()
            elif step_name == "phonon_dos":
                return self._run_phonon_dos()
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
        ensure_poscar(self.structure_file, self.relax_dir / "POSCAR")

        # 创建INCAR
        self._write_relax_incar(self.relax_dir / "INCAR")

        # 创建KPOINTS
        self._write_kpoints(self.relax_dir / "KPOINTS", self.relax_dir / "POSCAR", self.kspacing)

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
            poscar_relaxed = self.work_dir / "POSCAR_relaxed"
            shutil.copy(contcar, poscar_relaxed)
            self.steps_data['relaxed_structure'] = str(poscar_relaxed)

        logger.info("结构优化完成")
        return True

    def _prepare_phonon(self) -> bool:
        """Step 2: 准备声子计算"""
        logger.info("准备声子计算...")

        self.phonon_dir.mkdir(parents=True, exist_ok=True)

        # 使用优化后的结构
        relaxed_poscar = Path(self.steps_data.get('relaxed_structure', self.structure_file))

        # 复制到phonon目录并命名为POSCAR-init
        shutil.copy(relaxed_poscar, self.phonon_dir / "POSCAR-init")

        if self.method == "disp":
            # 使用phonopy生成位移超胞
            os.chdir(self.phonon_dir)

            supercell_str = f"{self.supercell[0]} {self.supercell[1]} {self.supercell[2]}"
            cmd = f"phonopy --dim=\"{supercell_str}\" -d -c POSCAR-init"

            logger.info(f"运行: {cmd}")
            ret = os.system(cmd)

            if ret != 0:
                logger.error("Phonopy生成位移超胞失败")
                return False

            # 检查生成的POSCAR-XXX文件
            poscar_files = list(self.phonon_dir.glob("POSCAR-[0-9]*"))
            logger.info(f"生成了 {len(poscar_files)} 个位移超胞")

            if len(poscar_files) == 0:
                logger.error("未生成位移超胞")
                return False

            self.steps_data['n_displacements'] = len(poscar_files)

        elif self.method == "dfpt":
            # DFPT法只需要一个超胞
            logger.info("使用DFPT方法，准备超胞...")
            # TODO: 使用phonopy或手动生成超胞

        logger.info("声子计算准备完成")
        return True

    def _run_phonon_calculate(self) -> bool:
        """Step 3: 执行声子计算"""
        logger.info("执行声子计算...")

        if self.method == "disp":
            # 位移法：需要计算所有disp-XXX
            n_disp = self.steps_data.get('n_displacements') or len(list(self.phonon_dir.glob("POSCAR-*")))
            if n_disp == 0:
                logger.error("未找到位移超胞，无法执行声子计算")
                return False

            # 1) 若未生成过输入，先准备所有位移输入
            if not self.steps_data.get("phonon_inputs_ready"):
                for i in range(1, n_disp + 1):
                    disp_num = str(i).zfill(3)
                    disp_dir = self.phonon_dir / f"disp-{disp_num}"
                    disp_dir.mkdir(exist_ok=True)

                    poscar_src = self.phonon_dir / f"POSCAR-{disp_num}"
                    shutil.copy(poscar_src, disp_dir / "POSCAR")

                    self._write_phonon_incar(disp_dir / "INCAR")
                    self._write_kpoints(disp_dir / "KPOINTS", disp_dir / "POSCAR", self.kspacing)

                    potcar_source = self.relax_dir / "POTCAR"
                    if potcar_source.exists():
                        shutil.copy(potcar_source, disp_dir / "POTCAR")
                    elif self.potcar_dir:
                        if not prepare_potcar(disp_dir / "POSCAR", self.potcar_dir, disp_dir / "POTCAR", self.potcar_type):
                            logger.error("POTCAR准备失败")
                            return False
                    else:
                        logger.warning(f"未找到POTCAR文件: {potcar_source}")

                    self._write_job_script(disp_dir, f"disp{disp_num}")

                self.steps_data["phonon_inputs_ready"] = True
                self._save_checkpoint()

            if self.prepare_only:
                logger.info("prepare_only=True，已生成所有位移输入和脚本，不提交。")
                return True

            # 2) 提交未完成的位移，已完成的跳过
            pending_jobs: List[tuple[str, Path]] = []
            for i in range(1, n_disp + 1):
                disp_num = str(i).zfill(3)
                disp_dir = self.phonon_dir / f"disp-{disp_num}"

                if check_vasp_completion(disp_dir):
                    logger.info(f"{disp_dir.name} 已完成，跳过提交")
                    continue

                job_script = disp_dir / f"run_disp{disp_num}.sh"
                job_id = self._submit_job(disp_dir, str(job_script))
                pending_jobs.append((job_id, disp_dir))
                logger.info(f"已提交位移计算 {i}/{n_disp}: {disp_dir.name}")

            # 3) 等待未完成的位移作业
            logger.info("等待所有未完成的位移计算完成...")
            for job_id, disp_dir in pending_jobs:
                if not self._wait_for_job(job_id, disp_dir, self.queue_system):
                    logger.error(f"位移计算失败: {disp_dir.name}")
                    return False

            logger.info("所有位移计算完成或已完成")

        elif self.method == "dfpt":
            # DFPT法
            # TODO: 实现DFPT计算
            logger.warning("DFPT方法尚未完全实现")

        return True

    def _run_phonon_band(self) -> bool:
        """Step 4: 后处理声子能带"""
        logger.info("后处理声子能带...")

        os.chdir(self.phonon_dir)

        if self.method == "disp":
            # 收集力常数
            n_disp = self.steps_data.get('n_displacements', 0)
            disp_str = "{" + f"001..{str(n_disp).zfill(3)}" + "}"

            cmd = f"phonopy -f disp-{disp_str}/vasprun.xml"
            logger.info(f"运行: {cmd}")
            os.system(cmd)

            # 生成能带配置
            self._write_phonon_band_conf(self.phonon_dir / "band.conf")

            # 计算声子能带
            cmd = "phonopy -p -s band.conf -c POSCAR-init"
            logger.info(f"运行: {cmd}")
            os.system(cmd)

            # 生成数据文件
            cmd = "phonopy-bandplot --gnuplot > band.dat"
            os.system(cmd)

        logger.info("声子能带后处理完成")
        return True

    def _run_phonon_dos(self) -> bool:
        """Step 5: 后处理声子DOS"""
        logger.info("后处理声子DOS...")

        os.chdir(self.phonon_dir)

        # 生成DOS配置
        self._write_phonon_dos_conf(self.phonon_dir / "mesh.conf")

        # 计算声子DOS
        cmd = "phonopy -p -t mesh.conf -c POSCAR-init"
        logger.info(f"运行: {cmd}")
        os.system(cmd)

        logger.info("声子DOS后处理完成")
        return True

    def _run_plotting(self) -> bool:
        """Step 6: 自动绘图"""
        logger.info("开始绘图...")

        self.plots_dir.mkdir(parents=True, exist_ok=True)

        try:
            # 绘制声子能带
            plotters.plot_phonon_band(
                self.phonon_dir,
                self.plots_dir / "phonon_band.png"
            )

            # 绘制声子DOS
            plotters.plot_phonon_dos(
                self.phonon_dir,
                self.plots_dir / "phonon_dos.png"
            )

            logger.info(f"所有图表已保存到: {self.plots_dir}")
            return True

        except Exception as e:
            logger.error(f"绘图失败: {e}", exc_info=True)
            return False

    def _write_relax_incar(self, incar_file: Path):
        """写入结构优化INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# Relaxation INCAR for Phonon\n")
            f.write("SYSTEM = Structure Relaxation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("EDIFF = 1E-8\n")  # 声子计算需要更高精度
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.01\n\n")
            f.write(f"PSTRESS = {self.pressure_kbar}\n\n")
            f.write("IBRION = 2\n")
            f.write("NSW = 200\n")
            f.write("ISIF = 3\n")
            f.write("EDIFFG = -1E-5\n\n")  # 更严格的收敛标准
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_phonon_incar(self, incar_file: Path):
        """写入声子计算INCAR"""
        with open(incar_file, 'w') as f:
            f.write("# Phonon INCAR (Displacement method)\n")
            f.write("SYSTEM = Phonon Calculation\n\n")
            f.write("PREC = Accurate\n")
            f.write("ENCUT = {}\n".format(self.encut if self.encut else 520))
            f.write("EDIFF = 1E-8\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.01\n")
            f.write("IBRION = -1\n")  # 单点能量计算
            f.write("NSW = 0\n")
            f.write(f"PSTRESS = {self.pressure_kbar}\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_kpoints(self, kpoints_file: Path, poscar_file: Path, kspacing: float):
        """写入K点（KSPACING -> Monkhorst-Pack 网格数）。"""
        from vasp.pipelines.utils import kspacing_to_mesh
        n1, n2, n3 = kspacing_to_mesh(poscar_file, kspacing)
        with open(kpoints_file, 'w') as f:
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{n1} {n2} {n3}\n")
            f.write("0 0 0\n")

    def _write_phonon_band_conf(self, conf_file: Path):
        """写入phonopy band.conf"""
        with open(conf_file, 'w') as f:
            f.write(f"DIM = {self.supercell[0]} {self.supercell[1]} {self.supercell[2]}\n")
            f.write("BAND = AUTO\n")
            f.write("BAND_POINTS = 101\n")

    def _write_phonon_dos_conf(self, conf_file: Path):
        """写入phonopy mesh.conf"""
        with open(conf_file, 'w') as f:
            f.write(f"DIM = {self.supercell[0]} {self.supercell[1]} {self.supercell[2]}\n")
            f.write("MP = 20 20 20\n")
            f.write("TPROP = T\n")

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
