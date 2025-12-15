import importlib.util
import logging
import os
import shlex
import subprocess
import time
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

_TOML_PATH = Path(__file__).resolve().parent.parent / "config" / "job_templates.toml"
_USER_TOML = Path.home() / ".config" / "vasp" / "job_templates.toml"
_LOCAL_TOML = Path.cwd() / "job_templates.local.toml"


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
) -> JobConfig:
    """
    读取 VASP 运行配置，优先级：工作目录 job_templates.local.toml > 用户级 ~/.config/vasp/job_templates.toml > 仓库默认 job_templates.toml。
    不再读取环境变量覆盖或 ~/.my_scriptrc.py 兼容。
    """
    cfg = {}

    # 1) 工作目录本地覆盖
    if _LOCAL_TOML.exists():
        defaults, templates = _load_templates_from_toml(_LOCAL_TOML)
    # 2) 用户级全局
    elif _USER_TOML.exists():
        defaults, templates = _load_templates_from_toml(_USER_TOML)
    # 3) 仓库内默认
    else:
        defaults, templates = _load_templates_from_toml(toml_path)

    if defaults.get("vasp_std"):
        cfg["vasp_std"] = defaults["vasp_std"]
    if defaults.get("vasp_gam"):
        cfg["vasp_gam"] = defaults["vasp_gam"]
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

    required = ["vasp_std", "vasp_gam"]
    missing = [k for k in required if not cfg.get(k)]
    if missing:
        raise ValueError(f"job_templates 配置缺少必需字段: {missing}，请在 job_templates.local.toml 或 ~/.config/vasp/job_templates.toml 中补充")

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
    workdir = script_path.parent
    jobid_file = workdir / "job_id.txt"

    # 如果已有作业ID且仍在运行，直接复用，不再重复提交
    if jobid_file.exists():
        try:
            existing_id = jobid_file.read_text().strip()
            if existing_id:
                from vasp.utils.job import is_job_active  # 避免循环导入
                active = is_job_active(existing_id, queue_system)
                if active:
                    logger.info(f"检测到目录已有运行中的任务 {existing_id}，跳过重复提交")
                    return existing_id
        except Exception:
            pass

    attempts = 0
    while True:
        attempts += 1
        try:
            if queue == "bash":
                subprocess.Popen(["bash", script_path.name], cwd=workdir)
                return "bash"
            if queue == "slurm":
                output = subprocess.check_output(
                    ["sbatch", script_path.name],
                    cwd=workdir,
                    text=True,
                    stderr=subprocess.STDOUT,
                ).strip()
                # 典型输出: "Submitted batch job 123456"
                parts = output.split()
                job_id = parts[-1] if parts else output
                jobid_file.write_text(str(job_id))
                return job_id
            if queue == "pbs":
                output = subprocess.check_output(["qsub", script_path.name], cwd=workdir, text=True).strip()
                jobid_file.write_text(str(output))
                return output
            if queue == "lsf":
                cmd = f"bsub < {shlex.quote(str(script_path))}"
                output = subprocess.check_output(["bash", "-lc", cmd], cwd=workdir, text=True).strip()
                # 典型输出: "Job <123456> is submitted ..."
                import re
                m = re.search(r"<(\d+)>", output)
                job_id = m.group(1) if m else output
                jobid_file.write_text(str(job_id))
                return job_id
        except subprocess.CalledProcessError as exc:
            msg = exc.output if hasattr(exc, "output") else str(exc)
            if queue == "slurm" and "AssocMaxSubmitJobLimit" in msg:
                logger.warning("提交超出作业数限制，等待20秒后重试（第 %d 次）...", attempts)
                time.sleep(20)
                continue
            logger.error("队列提交失败，请检查集群环境或脚本", exc_info=exc)
            raise
        except Exception as exc:  # pragma: no cover - 运行时容错
            logger.error("队列提交失败，请检查集群环境或脚本", exc_info=exc)
            raise

    # 未知队列，回退（理论不会到这里）
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
