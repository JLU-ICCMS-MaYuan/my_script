import importlib.util
import logging
import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

try:  # Python 3.11+
    import tomllib  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - 兼容旧版
    tomllib = None

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

_TOML_PATH = Path(__file__).resolve().parent.parent / "config" / "job_templates.toml"


def _load_templates_from_toml(toml_path: Path) -> tuple[dict, dict]:
    """从 TOML 配置读取 defaults 与 queue header 模板。"""
    if not tomllib or not toml_path.exists():
        return {}, {}

    try:
        with toml_path.open("rb") as f:
            data = tomllib.load(f)
        defaults = data.get("defaults", {}) or {}
        templates_raw = data.get("templates", {}) or {}
        templates = {k.lower(): (v.get("header") or "").strip() for k, v in templates_raw.items()}
        return defaults, templates
    except Exception as exc:  # pragma: no cover
        logger.warning("读取 job_templates.toml 失败，使用默认模板", exc_info=exc)
        return {}, {}


def load_job_config(
    toml_path: Path = _TOML_PATH,
    rc_path: Path = Path.home() / ".my_scriptrc.py",
) -> JobConfig:
    """
    读取 VASP 运行配置，优先级：环境变量 > TOML 模板 > ~/.my_scriptrc.py（兼容）> 内置默认。

    环境变量：
    - VASP_STD / VASP_GAM / POTCAR_DIR / VASP_MPI_PROCS
    - JOB_HEADER_BASH / JOB_HEADER_SLURM / JOB_HEADER_PBS / JOB_HEADER_LSF
    """
    cfg = _DEFAULTS.copy()

    # 兼容旧 rc（低优先级）
    rc_file = Path(rc_path)
    if rc_file.exists():
        try:
            spec = importlib.util.spec_from_file_location("_user_rc", rc_file)
            module = importlib.util.module_from_spec(spec)
            assert spec and spec.loader
            spec.loader.exec_module(module)  # type: ignore[attr-defined]
            rc_values = {
                "vasp_std": getattr(module, "vaspstd_path", None),
                "vasp_gam": getattr(module, "vaspgam_path", None),
                "potcar_dir": getattr(module, "potcar_dir", None),
                "bashtitle": getattr(module, "bashtitle", None),
                "slurmtitle": getattr(module, "slurmtitle", None),
                "pbstitle": getattr(module, "pbstitle", None),
                "lsftitle": getattr(module, "lsftitle", None),
                "default_mpi_procs": getattr(module, "default_mpi_procs", None),
            }
            for k, v in rc_values.items():
                if v:
                    cfg[k] = v
        except Exception as exc:  # pragma: no cover
            logger.warning("加载 ~/.my_scriptrc.py 失败，跳过 rc", exc_info=exc)

    # TOML 模板（中优先级）
    defaults, templates = _load_templates_from_toml(toml_path)
    if defaults.get("potcar_dir"):
        cfg["potcar_dir"] = defaults["potcar_dir"]
    if defaults.get("mpi_procs"):
        cfg["default_mpi_procs"] = defaults["mpi_procs"]

    for queue, header in templates.items():
        if queue == "bash":
            cfg["bashtitle"] = header or cfg["bashtitle"]
        if queue == "slurm":
            cfg["slurmtitle"] = header or cfg["slurmtitle"]
        if queue == "pbs":
            cfg["pbstitle"] = header or cfg["pbstitle"]
        if queue == "lsf":
            cfg["lsftitle"] = header or cfg["lsftitle"]

    # 环境变量（最高优先级）
    env_overrides = {
        "vasp_std": os.environ.get("VASP_STD"),
        "vasp_gam": os.environ.get("VASP_GAM"),
        "potcar_dir": os.environ.get("POTCAR_DIR"),
        "default_mpi_procs": os.environ.get("VASP_MPI_PROCS"),
        "bashtitle": os.environ.get("JOB_HEADER_BASH"),
        "slurmtitle": os.environ.get("JOB_HEADER_SLURM"),
        "pbstitle": os.environ.get("JOB_HEADER_PBS"),
        "lsftitle": os.environ.get("JOB_HEADER_LSF"),
    }
    for k, v in env_overrides.items():
        if v:
            cfg[k] = v

    potcar_dir = Path(cfg["potcar_dir"]) if cfg.get("potcar_dir") else None

    return JobConfig(
        vasp_std=str(cfg["vasp_std"]),
        vasp_gam=str(cfg["vasp_gam"]),
        potcar_dir=potcar_dir,
        bashtitle=str(cfg["bashtitle"]).strip(),
        slurmtitle=str(cfg["slurmtitle"]).strip() if cfg.get("slurmtitle") else None,
        pbstitle=str(cfg["pbstitle"]).strip() if cfg.get("pbstitle") else None,
        lsftitle=str(cfg["lsftitle"]).strip() if cfg.get("lsftitle") else None,
        default_mpi_procs=int(cfg.get("default_mpi_procs") or 8),
    )


def select_job_header(queue_system: str, cfg: JobConfig) -> str:
    """脚本头选择器，所有模块共用，按队列类型选择 header。"""
    header_map = {
        "bash": cfg.bashtitle,
        "slurm": cfg.slurmtitle or cfg.bashtitle,
        "pbs": cfg.pbstitle or cfg.bashtitle,
        "lsf": cfg.lsftitle or cfg.bashtitle,
    }
    header = header_map.get((queue_system or "bash"), cfg.bashtitle) or cfg.bashtitle
    return header.strip()


def write_job_script(
    work_dir: Path,
    job_name: str,
    queue_system: str,
    cfg: JobConfig,
    use_gamma: bool = False,
    mpi_procs: Optional[Union[int, str]] = None,
) -> Path:
    """根据 queue_system 生成作业脚本。"""
    work_dir = Path(work_dir)
    queue = (queue_system or "bash").lower()
    script_file = work_dir / f"run_{job_name}.sh"

    header = select_job_header(queue, cfg)
    binary = cfg.vasp_gam if use_gamma else cfg.vasp_std
    # 支持字符串或数字：字符串直接作为启动命令前缀，数字走 mpirun -np
    run_line: str
    if mpi_procs is None:
        mpi = cfg.default_mpi_procs or 8
        run_line = f"mpirun -np {mpi} {binary} > vasp.log"
    elif isinstance(mpi_procs, str):
        mpi_str = mpi_procs.strip()
        if mpi_str.isdigit():
            run_line = f"mpirun -np {mpi_str} {binary} > vasp.log"
        else:
            run_line = f"{mpi_str} {binary} > vasp.log"
    else:
        run_line = f"mpirun -np {mpi_procs} {binary} > vasp.log"

    lines = [header, run_line]

    script_file.write_text("\n".join(lines) + "\n")
    script_file.chmod(0o755)
    return script_file


def submit_job(script_path: Path, queue_system: str) -> str:
    """提交作业，不同队列系统自动选择命令；失败时回落到本地 bash 执行。"""
    queue = (queue_system or "bash").lower()
    script_path = Path(script_path).resolve()

    try:
        if queue == "bash":
            subprocess.Popen(["bash", script_path.name], cwd=script_path.parent)
            return "bash"
        if queue == "slurm":
            output = subprocess.check_output(["sbatch", str(script_path)], text=True).strip()
            # 典型输出: "Submitted batch job 123456"
            parts = output.split()
            return parts[-1] if parts else output
        if queue == "pbs":
            output = subprocess.check_output(["qsub", str(script_path)], text=True).strip()
            return output
        if queue == "lsf":
            cmd = f"bsub < {shlex.quote(str(script_path))}"
            output = subprocess.check_output(["bash", "-lc", cmd], text=True).strip()
            # 典型输出: "Job <123456> is submitted ..."
            import re
            m = re.search(r"<(\d+)>", output)
            return m.group(1) if m else output
    except Exception as exc:  # pragma: no cover - 运行时容错
        logger.warning("队列提交失败，改为本地 bash 运行", exc_info=exc)
        subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
        return "bash_fallback"

    # 未知队列，回退
    logger.warning("未知队列系统 %s，回退到 bash", queue_system)
    subprocess.Popen(["bash", str(script_path)], cwd=script_path.parent)
    return "bash_unknown"


def is_job_active(job_id: str, queue_system: Optional[str]) -> Optional[bool]:
    """
    查询作业是否仍在队列中。

    返回值：True=仍在运行/排队，False=不在队列，None=无法判断或使用bash。
    """
    queue = (queue_system or "bash").lower()

    if queue == "bash":
        return None

    try:
        if queue == "slurm":
            output = subprocess.check_output(["squeue", "-j", str(job_id)], text=True)
            # squeue 输出含表头 + 任务行
            return len(output.strip().splitlines()) > 1

        if queue == "pbs":
            output = subprocess.check_output(["qstat", str(job_id)], text=True)
            return job_id in output

        if queue == "lsf":
            output = subprocess.check_output(["bjobs", str(job_id)], text=True)
            return job_id in output
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        logger.warning("未找到队列命令，跳过队列状态检查")
        return None
    except Exception as exc:  # pragma: no cover
        logger.warning("队列状态检查异常，跳过", exc_info=exc)
        return None

    return None
