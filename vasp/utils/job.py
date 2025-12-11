import importlib.util
import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class JobConfig:
    """封装作业提交相关的用户配置。"""

    vasp_std: str
    vasp_gam: str
    potcar_dir: Optional[Path]
    bashtitle: str
    slurmtitle: Optional[str]
    pbstitle: Optional[str]
    lsftitle: Optional[str]
    default_mpi_procs: int = 8


_DEFAULTS = {
    "vasp_std": "/home/mayuan/soft/vasp.6.3.2/bin/vasp_std",
    "vasp_gam": "/home/mayuan/soft/vasp.6.3.2/bin/vasp_gam",
    "potcar_dir": "/home/mayuan/pot/vasp_pot/potpaw_PBE54",
    "bashtitle": "#!/bin/sh\nsource /home/mayuan/intel/oneapi/setvars.sh\nulimit -s unlimited\n",
    "slurmtitle": "#!/bin/sh\n#SBATCH  --job-name=vasp\n#SBATCH  --output=log.out\n#SBATCH  --error=log.err\n##SBATCH  --partition=intel6240r_384\n#SBATCH  --partition=intel6240r_192\n#SBATCH  --nodes=1\n#SBATCH  --ntasks=48\n#SBATCH  --ntasks-per-node=48\n#SBATCH  --cpus-per-task=1\n\nsource /home/mayuan/intel/oneapi/setvars.sh\nulimit -s unlimited\nexport I_MPI_ADJUST_REDUCE=3\nexport MPIR_CVAR_COLL_ALIAS_CHECK=0\n",
    "pbstitle": "#!/bin/sh\n#PBS -N    mayqe\n#PBS -q    liuhy\n#PBS -l    nodes=1:ppn=28\n#PBS -j    oe\n#PBS -V\n\nsource /home/mayuan/intel/oneapi/setvars.sh\nulimit -s unlimited\ncd $PBS_O_WORKDIR\n",
    "lsftitle": "#!/bin/bash\n#BSUB -n 56\n#BSUB -q normal\n#BSUB -J myjob\n#BSUB -R 'span[ptile=56]'\n#BSUB -o operation.log\n\nsource /data/env/inteloneapi2021\nulimit -s unlimited\n",
    "default_mpi_procs": 8,
}


def load_job_config(rc_path: Path = Path.home() / ".my_scriptrc.py") -> JobConfig:
    """
    加载 ~/.my_scriptrc.py 中的 VASP 配置；缺失时使用默认值。

    用户可在 rc 中提供：
    - vaspstd_path, vaspgam_path, potcar_dir
    - bashtitle, slurmtitle, pbstitle, lsftitle
    """
    cfg = _DEFAULTS.copy()
    rc_file = Path(rc_path)

    if rc_file.exists():
        try:
            spec = importlib.util.spec_from_file_location("_user_rc", rc_file)
            module = importlib.util.module_from_spec(spec)
            assert spec and spec.loader
            spec.loader.exec_module(module)  # type: ignore[attr-defined]

            cfg["vasp_std"] = getattr(module, "vaspstd_path", cfg["vasp_std"])
            cfg["vasp_gam"] = getattr(module, "vaspgam_path", cfg["vasp_gam"])
            cfg["potcar_dir"] = getattr(module, "potcar_dir", cfg["potcar_dir"])
            cfg["bashtitle"] = getattr(module, "bashtitle", cfg["bashtitle"])
            cfg["slurmtitle"] = getattr(module, "slurmtitle", cfg["slurmtitle"])
            cfg["pbstitle"] = getattr(module, "pbstitle", cfg["pbstitle"])
            cfg["lsftitle"] = getattr(module, "lsftitle", cfg["lsftitle"])
            cfg["default_mpi_procs"] = getattr(module, "default_mpi_procs", cfg["default_mpi_procs"])
        except Exception as exc:  # pragma: no cover - 容错为主
            logger.warning("加载 ~/.my_scriptrc.py 失败，使用默认配置", exc_info=exc)
    else:
        logger.info("未找到 ~/.my_scriptrc.py，使用默认 VASP 配置")

    potcar_dir = Path(cfg["potcar_dir"]) if cfg.get("potcar_dir") else None

    return JobConfig(
        vasp_std=str(cfg["vasp_std"]),
        vasp_gam=str(cfg["vasp_gam"]),
        potcar_dir=potcar_dir,
        bashtitle=str(cfg["bashtitle"]).strip(),
        slurmtitle=str(cfg["slurmtitle"]).strip() if cfg.get("slurmtitle") else None,
        pbstitle=str(cfg["pbstitle"]).strip() if cfg.get("pbstitle") else None,
        lsftitle=str(cfg["lsftitle"]).strip() if cfg.get("lsftitle") else None,
        default_mpi_procs=int(cfg.get("default_mpi_procs", 8)),
    )


def _select_header(queue_system: str, cfg: JobConfig) -> str:
    header_map = {
        "bash": cfg.bashtitle,
        "slurm": cfg.slurmtitle or cfg.bashtitle,
        "pbs": cfg.pbstitle or cfg.bashtitle,
        "lsf": cfg.lsftitle or cfg.bashtitle,
    }
    header = header_map.get(queue_system, cfg.bashtitle) or cfg.bashtitle
    return header.strip()


def write_job_script(
    work_dir: Path,
    job_name: str,
    queue_system: str,
    cfg: JobConfig,
    use_gamma: bool = False,
    mpi_procs: Optional[int] = None,
) -> Path:
    """根据 queue_system 生成作业脚本。"""
    work_dir = Path(work_dir)
    queue = (queue_system or "bash").lower()
    script_file = work_dir / f"run_{job_name}.sh"

    header = _select_header(queue, cfg)
    binary = cfg.vasp_gam if use_gamma else cfg.vasp_std
    mpi = mpi_procs or cfg.default_mpi_procs or 8

    lines = [
        header,
        f"cd {work_dir}",
        f"mpirun -np {mpi} {binary} > vasp.log",
    ]

    script_file.write_text("\n".join(lines) + "\n")
    script_file.chmod(0o755)
    return script_file


def submit_job(script_path: Path, queue_system: str) -> str:
    """提交作业，不同队列系统自动选择命令；失败时回落到本地 bash 执行。"""
    queue = (queue_system or "bash").lower()
    script_path = Path(script_path)

    try:
        if queue == "bash":
            subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
            return "bash"
        if queue == "slurm":
            output = subprocess.check_output(["sbatch", str(script_path)], text=True)
            return output.strip()
        if queue == "pbs":
            output = subprocess.check_output(["qsub", str(script_path)], text=True)
            return output.strip()
        if queue == "lsf":
            cmd = f"bsub < {shlex.quote(str(script_path))}"
            output = subprocess.check_output(["bash", "-lc", cmd], text=True)
            return output.strip()
    except Exception as exc:  # pragma: no cover - 运行时容错
        logger.warning("队列提交失败，改为本地 bash 运行", exc_info=exc)
        subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
        return "bash_fallback"

    # 未知队列，回退
    logger.warning("未知队列系统 %s，回退到 bash", queue_system)
    subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
    return "bash_unknown"
