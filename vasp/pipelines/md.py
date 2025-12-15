import logging
import shutil
from pathlib import Path
from typing import Optional, List

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import prepare_potcar, ensure_poscar
from vasp.utils.job import load_job_config, write_job_script, submit_job

logger = logging.getLogger(__name__)


class MdPipeline(BasePipeline):
    """分子动力学（NVT）Pipeline。"""

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        potim: float = 1.0,
        tebeg: float = 300.0,
        teend: float = 300.0,
        nsw: int = 200,
        kspacing: float = 0.2,
        encut: Optional[float] = None,
        queue_system: Optional[str] = None,
        mpi_procs: Optional[str] = None,
        potcar_dir: Optional[Path] = None,
        potcar_type: str = "PBE",
        include_relax: bool = True,
        custom_steps: Optional[List[str]] = None,
        pressure: float = 0.0,
        **kwargs,
    ):
        super().__init__(structure_file, work_dir, pressure=pressure, **kwargs)

        self.job_cfg = load_job_config()
        self.potim = potim
        self.tebeg = tebeg
        self.teend = teend
        self.nsw = nsw
        self.kspacing = kspacing
        self.encut = encut
        self.queue_system = queue_system or "bash"
        self.mpi_procs = mpi_procs
        default_potcar = self.job_cfg.potcar_dir if self.job_cfg else None
        self.potcar_dir = Path(potcar_dir) if potcar_dir else default_potcar
        self.potcar_type = potcar_type
        self.include_relax = include_relax
        self.custom_steps = self._normalize_steps(custom_steps)
        self.pressure = pressure

        self.relax_dir = self.work_dir / "01_relax" if self.include_relax else None
        self.md_dir = self.work_dir / ("02_md" if self.include_relax else "01_md")

    def get_steps(self):
        if self.custom_steps:
            return self.custom_steps
        return ["relax", "md"]

    def execute_step(self, step_name: str) -> bool:
        if step_name == "relax":
            return self._run_relax()
        if step_name == "md":
            return self._run_md()
        logger.error(f"未知步骤: {step_name}")
        return False

    def _normalize_steps(self, custom_steps: Optional[List[str]]) -> Optional[List[str]]:
        if not custom_steps:
            return None
        allowed = ["relax", "md"]
        normalized: List[str] = []
        for step in custom_steps:
            name = step.strip().lower()
            if name in allowed:
                if name not in normalized:
                    normalized.append(name)
            else:
                logger.warning(f"忽略未支持的步骤: {name}")
        if self.include_relax and "md" in normalized and "relax" not in normalized:
            normalized.insert(0, "relax")
        return normalized or None

    def _run_relax(self) -> bool:
        logger.info("执行结构优化(为MD准备)...")

        if not self.relax_dir:
            logger.error("未启用 include_relax，不能运行 relax 步骤")
            return False

        self.relax_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(self.structure_file, self.relax_dir / "POSCAR")

        self._write_relax_incar(self.relax_dir / "INCAR")
        self._write_kpoints(self.relax_dir / "KPOINTS", self.relax_dir / "POSCAR", self.kspacing)

        if self.potcar_dir:
            if not prepare_potcar(
                self.relax_dir / "POSCAR",
                self.potcar_dir,
                self.relax_dir / "POTCAR",
                self.potcar_type,
            ):
                logger.error("POTCAR准备失败")
                return False
        else:
            logger.warning("未提供 potcar_dir，请确保已手动准备 POTCAR")

        job_script = self._write_job_script(self.relax_dir, "relax")

        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.relax_dir, job_script)

        if not self._wait_for_job(job_id, self.relax_dir, self.queue_system):
            return False

        if not self._check_convergence(self.relax_dir):
            logger.error("结构优化未收敛")
            return False

        contcar = self.relax_dir / "CONTCAR"
        if contcar.exists():
            relaxed = self.work_dir / "POSCAR_relaxed"
            shutil.copy(contcar, relaxed)
            self.steps_data["relaxed_structure"] = str(relaxed)

        logger.info("结构优化完成")
        return True

    def _run_md(self) -> bool:
        logger.info("执行分子动力学计算...")

        self.md_dir.mkdir(parents=True, exist_ok=True)
        # 使用relax结果优先
        source_poscar = Path(self.structure_file)
        if self.include_relax and self.steps_data.get("relaxed_structure"):
            source_poscar = Path(self.steps_data["relaxed_structure"])
        ensure_poscar(source_poscar, self.md_dir / "POSCAR")

        self._write_md_incar(self.md_dir / "INCAR")
        self._write_kpoints(self.md_dir / "KPOINTS", self.md_dir / "POSCAR", self.kspacing)

        if self.include_relax and self.relax_dir and (self.relax_dir / "POTCAR").exists():
            shutil.copy(self.relax_dir / "POTCAR", self.md_dir / "POTCAR")
        elif self.potcar_dir:
            if not prepare_potcar(
                self.md_dir / "POSCAR",
                self.potcar_dir,
                self.md_dir / "POTCAR",
                self.potcar_type,
            ):
                logger.error("POTCAR准备失败")
                return False
        else:
            logger.warning("未提供 potcar_dir，请确保已手动准备 POTCAR")

        job_script = self._write_job_script(self.md_dir, "md")

        if self.prepare_only:
            logger.info("prepare_only=True，仅生成输入和脚本，不提交。")
            return True

        job_id = self._submit_job(self.md_dir, job_script)

        if not self._wait_for_job(job_id, self.md_dir, self.queue_system):
            return False

        if not self._check_job_completed(self.md_dir):
            logger.error("MD 计算未完成")
            return False

        logger.info("分子动力学计算完成")
        return True

    def _write_md_incar(self, incar_file: Path):
        """写入 MD INCAR。"""
        with open(incar_file, "w") as f:
            f.write("# Molecular Dynamics INCAR\n")
            f.write("SYSTEM = MD Simulation\n\n")
            f.write("PREC = Accurate\n")
            f.write(f"ENCUT = {self.encut if self.encut else 520}\n")
            f.write("EDIFF = 1E-6\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n\n")
            f.write("IBRION = 0\n")
            f.write("MDALGO = 2\n")
            f.write(f"NSW = {self.nsw}\n")
            f.write(f"POTIM = {self.potim}\n")
            f.write(f"TEBEG = {self.tebeg}\n")
            f.write(f"TEEND = {self.teend}\n")
            f.write("SMASS = 0\n")
            f.write(f"PSTRESS = {self.pressure_kbar}\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_relax_incar(self, incar_file: Path):
        """写入用于MD前结构优化的 INCAR。"""
        with open(incar_file, "w") as f:
            f.write("# Relaxation for MD\n")
            f.write("SYSTEM = Relax before MD\n\n")
            f.write("PREC = Accurate\n")
            f.write(f"ENCUT = {self.encut if self.encut else 520}\n")
            f.write("EDIFF = 1E-6\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n")
            f.write("IBRION = 2\n")
            f.write("NSW = 120\n")
            f.write("ISIF = 3\n")
            f.write("EDIFFG = -0.02\n")
            f.write(f"PSTRESS = {self.pressure_kbar}\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_kpoints(self, kpoints_file: Path, poscar_file: Path, kspacing: float):
        from vasp.pipelines.utils import kspacing_to_mesh
        n1, n2, n3 = kspacing_to_mesh(poscar_file, kspacing)
        with open(kpoints_file, "w") as f:
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{n1} {n2} {n3}\n")
            f.write("0 0 0\n")

    def _write_job_script(self, work_dir: Path, job_name: str) -> str:
        script_path = write_job_script(
            work_dir=work_dir,
            job_name=job_name,
            queue_system=self.queue_system,
            cfg=self.job_cfg,
            mpi_procs=self.mpi_procs,
        )
        return str(script_path)

    def _submit_job(self, work_dir: Path, job_script: str) -> str:
        return submit_job(Path(job_script), self.queue_system)
