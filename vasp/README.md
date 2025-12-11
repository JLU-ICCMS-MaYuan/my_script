# VASP 计算工具 - 使用与开发指引

## 配置来源
- 默认读取 `~/.my_scriptrc.py` 中的 `vaspstd_path/vaspgam_path/potcar_dir` 以及 `bashtitle/slurmtitle/pbstitle/lsftitle` 作为作业脚本头；缺失时回退到 `vasp/utils/job.py` 内置路径（同样基于 rc 提供的示例）。
- 队列系统通过 `-j/--job-system` 选择 `bash/slurm/pbs/lsf`；MPI 进程数可用 `--mpi-procs` 指定（默认 8 或 rc 中的设置）。
- POTCAR 默认读取 rc 的 `potcar_dir`，也可用 `--potcar-dir` 覆盖；`--potcar-type` 选择 PBE/LDA/PW91。

## 功能与流程
- `relax`：仅结构优化（步骤目录 `01_relax`）。可批量、并行。
- `electronic`：自动串联 `01_relax → 02_scf → 03_dos → 04_band → 05_elf → 06_cohp → plots`，并生成图像。
- `phonon`：支持 `--method disp/dfpt`；位移法自动生成超胞并批量提交 `disp-***` 任务，后处理能带/DOS并绘图。
- `md`：NVT 分子动力学（`MDALGO=2`），仅单任务流程。
- 所有 Pipeline 支持断点续跑与状态记录：`pipeline_checkpoint.json`、`pipeline_report.txt`。

## 快速示例
```bash
# 结构优化（使用 rc 中的 potcar/脚本头）
vasp relax -i POSCAR -w ./relax -j slurm --mpi-procs 48

# 电子性质全流程（含 ELF/COHP）
vasp electronic -i POSCAR -w ./electronic \
  --kspacing 0.2 --encut 500 --include-elf --include-cohp \
  --potcar-type PBE -j pbs

# 声子（位移法，超胞 2×2×2）
vasp phonon -i POSCAR -w ./phonon --supercell 2 2 2 \
  --method disp --kspacing 0.3 -j bash

# 分子动力学
vasp md -i POSCAR -w ./md --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf

# 批量并行（自动扫描结构目录）
vasp electronic -i ./structures -w ./batch_elec --parallel --max-workers 4
```

## 批量与并行
- 目录输入或 `--batch` 触发批量模式；`--parallel --max-workers N` 控制并行。
- 每个结构生成独立子目录并统计成功/失败，汇总写入 `batch_summary.txt`。

## 工作目录与产物
- 电子流程：分步目录与 `plots/`（band/dos/elf/cohp）；
- 声子流程：`02_phonon/disp-*` + phonopy 后处理产物与 `plots/`；
- 断点文件 `pipeline_checkpoint.json` 支持续跑，报告 `pipeline_report.txt` 便于回溯。

## 注意事项
- 队列提交通过 rc/默认脚本头生成并调用 `mpirun`；命令缺失时回退本地 bash 执行，bash 环境即可直接跑。
- `_wait_for_job` 同时检查队列状态（slurm/pbs/lsf）与 OUTCAR 关键字，若队列中不存在该任务会提前报错；可根据集群行为调整判据或超时。
- 仓库不内置 VASP 可执行文件/赝势，请确保路径、权限与 MPI 环境正确。
