import logging
import shutil
from pathlib import Path
from typing import Optional

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import prepare_potcar
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
        mpi_procs: Optional[int] = None,
        potcar_dir: Optional[Path] = None,
        potcar_type: str = "PBE",
        **kwargs,
    ):
        super().__init__(structure_file, work_dir, **kwargs)

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

        self.md_dir = self.work_dir / "01_md"

    def get_steps(self):
        return ["md"]

    def execute_step(self, step_name: str) -> bool:
        if step_name == "md":
            return self._run_md()
        logger.error(f"未知步骤: {step_name}")
        return False

    def _run_md(self) -> bool:
        logger.info("执行分子动力学计算...")

        self.md_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(self.structure_file, self.md_dir / "POSCAR")

        self._write_md_incar(self.md_dir / "INCAR")
        self._write_kpoints(self.md_dir / "KPOINTS", self.kspacing)

        if self.potcar_dir:
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
        job_id = self._submit_job(self.md_dir, job_script)

        if self.submit_only:
            logger.info("submit_only=True，已提交 md 作业，退出等待。")
            return True

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
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .FALSE.\n")

    def _write_kpoints(self, kpoints_file: Path, kspacing: float):
        with open(kpoints_file, "w") as f:
            f.write("Automatic mesh\n")
            f.write("0\n")
            f.write("Gamma\n")
            f.write(f"{kspacing} {kspacing} {kspacing}\n")

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
