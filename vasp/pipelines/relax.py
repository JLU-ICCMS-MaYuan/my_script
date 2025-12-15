import logging
import shutil
from pathlib import Path
from typing import Optional

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import prepare_potcar, ensure_poscar, find_symmetry
from vasp.utils.job import load_job_config, write_job_script, submit_job

logger = logging.getLogger(__name__)


class RelaxPipeline(BasePipeline):
    """仅执行结构优化的 Pipeline。"""

    def __init__(
        self,
        structure_file: Path,
        work_dir: Path,
        kspacing: float = 0.2,
        encut: Optional[float] = None,
        queue_system: Optional[str] = None,
        mpi_procs: Optional[str] = None,
        potcar_dir: Optional[Path] = None,
        potcar_type: str = "PBE",
        pressure: float = 0.0,
        **kwargs,
    ):
        super().__init__(structure_file, work_dir, pressure=pressure, **kwargs)

        self.job_cfg = load_job_config()
        self.kspacing = kspacing
        self.encut = encut
        self.queue_system = queue_system or "bash"
        self.mpi_procs = mpi_procs
        default_potcar = self.job_cfg.potcar_dir if self.job_cfg else None
        self.potcar_dir = Path(potcar_dir) if potcar_dir else default_potcar
        self.potcar_type = potcar_type
        self.pressure = pressure

        self.relax_dir = self.work_dir / "01_relax"

    def get_steps(self):
        return ["relax"]

    def execute_step(self, step_name: str) -> bool:
        if step_name == "relax":
            return self._run_relax()
        logger.error(f"未知步骤: {step_name}")
        return False

    def _run_relax(self) -> bool:
        logger.info("执行结构优化...")

        self.relax_dir.mkdir(parents=True, exist_ok=True)
        ensure_poscar(self.structure_file, self.relax_dir / "POSCAR")

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

            # 对称性分析，生成原胞/标准晶胞（放在压强目录）
            prim, std, sg = find_symmetry(relaxed, self.work_dir, symprec=1e-3)
            if prim:
                self.steps_data["primitive_structure"] = str(prim)
            if std:
                self.steps_data["conventional_structure"] = str(std)
            if sg:
                self.steps_data["spacegroup"] = sg
            self._save_checkpoint()

        logger.info("结构优化完成")
        return True

    def _write_relax_incar(self, incar_file: Path):
        """写入结构优化 INCAR。"""
        with open(incar_file, "w") as f:
            f.write("# Relaxation INCAR\n")
            f.write("SYSTEM = Structure Relaxation\n\n")
            f.write("PREC = Accurate\n")
            f.write(f"ENCUT = {self.encut if self.encut else 520}\n")
            f.write("EDIFF = 1E-6\n")
            f.write("ISMEAR = 0\n")
            f.write("SIGMA = 0.05\n\n")
            f.write(f"PSTRESS = {self.pressure_kbar}\n\n")
            f.write("IBRION = 2\n")
            f.write("NSW = 200\n")
            f.write("ISIF = 3\n")
            f.write("EDIFFG = -0.01\n\n")
            f.write("LWAVE = .FALSE.\n")
            f.write("LCHARG = .TRUE.\n")

    def _write_kpoints(self, kpoints_file: Path, poscar_file: Path, kspacing: float):
        """写入自动 K 点（由 KSPACING 计算网格数）。"""
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
