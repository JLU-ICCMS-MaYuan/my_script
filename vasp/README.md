# VASP 计算工具 - 使用与开发指引

## 配置来源
- 程序路径与赝势：优先读取环境变量 `VASP_STD`、`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`；其次读取 `vasp/config/job_templates.toml` 的 defaults；最后回退到内置默认（兼容旧 `.my_scriptrc.py`）。
- 队列脚本头：`vasp/config/job_templates.toml` 中按 bash/slurm/pbs/lsf 提供模板，可用环境变量 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖。
- 队列选择：`-j/--job-system` 取 `bash/slurm/pbs/lsf`；MPI 进程数可用 `--mpi-procs` 指定（默认 8 或配置值）。
- POTCAR：默认使用 `POTCAR_DIR` 或模板 defaults，可被 `--potcar-dir` 覆盖；`--potcar-type` 选择 PBE/LDA/PW91。

## 执行模式（统一约定）
- 默认仅准备：不加 `--submit` 时只生成当前待执行步骤的输入与 `run_*.sh`，不提交、不等待；多步流程需多次调用推进下一步（断点文件自动跳过已完成步骤）。
- 提交但不等待：加 `--submit` 时提交“当前待办步骤”后立即返回。
- Master 脚本（电子流程）：`--master-script [--steps ...]` 生成 `run_master.sh`，按所选步骤依次执行并自动拷贝依赖；再加 `--submit` 只提交这个脚本。
- Bash/队列：`-j bash/slurm/pbs/lsf` 选择脚本头，`--mpi-procs` 控制 `mpirun -np`。

## 各命令的执行流程
- 结构优化 `vasp relax`：生成 `01_relax`（POSCAR/INCAR/KPOINTS/POTCAR + run_relax.sh）→ 提交/准备。
- 电子性质 `vasp electronic`：
  - 默认步骤 `01_relax → 02_scf → 03_dos → 04_band → [05_elf] → [06_cohp] → plots`。
  - 每次调用只处理当前待办一步；`--steps` 可定制子集（如 `relax,scf,elf`）；`--master-script` 生成一键串行脚本。
- 声子 `vasp phonon`：
  - 默认仅声子：`02_phonon` 下生成位移超胞（disp）或 DFPT 输入 → 对每个 `disp-*` 准备/提交 → phonopy 后处理能带/DOS → 绘图。
  - 若需先优化：用 `vasp relax phonon ...` 或 `vasp phonon --with-relax ...`，流程变为 `01_relax` 后再进入声子步骤。
- 分子动力学 `vasp md`：
  - 默认仅 MD：生成 `01_md` 输入与脚本。
  - 若需先优化：`vasp relax md ...` 或 `vasp md --with-relax ...`，先 `01_relax`，再 `02_md`，并自动继承 CONTCAR/POTCAR。
- 批量模式：输入为目录或显式 `--batch` 时，扫描结构文件，逐个创建子目录独立运行；可配合 `--parallel --max-workers N` 并行。

## 快速示例与步骤说明
- 结构优化（1 步：`01_relax`，准备 INCAR/KPOINTS/POTCAR + `run_relax.sh`）
  ```bash
  vasp relax -i POSCAR -w ./relax -j slurm --mpi-procs 48          # 仅准备
  vasp relax -i POSCAR -w ./relax -j slurm --mpi-procs 48 --submit # 提交 relax
  ```
- 电子性质（6~7 步：`01_relax → 02_scf → 03_dos → 04_band → [05_elf] → [06_cohp] → plots`）
  ```bash
  vasp electronic -i POSCAR -w ./electronic --kspacing 0.2 --encut 500 \
    --include-elf --include-cohp --potcar-type PBE -j pbs                # 准备当前步
  vasp electronic -i POSCAR -w ./electronic --kspacing 0.2 --encut 500 \
    --include-elf --include-cohp --potcar-type PBE -j pbs --submit       # 提交当前步
  # 提示：每次运行只处理一个待办步骤；当 relax 完成后再运行同一命令推进到 scf，依次推进到 band/elf/cohp，最后自动绘图。
  # 一次提交串行跑完（例如 relax→scf→elf），并生成 plots：
  vasp electronic -i POSCAR -w ./electronic --kspacing 0.2 --steps relax,scf,elf \
    --master-script -j pbs                       # 生成 run_master.sh（不提交）
  vasp electronic -i POSCAR -w ./electronic --kspacing 0.2 --steps relax,scf,elf \
    --master-script -j pbs --submit              # 提交 run_master.sh，脚本内依次执行各步
  ```
- 声子（位移法典型 3 段：`01_relax → 02_phonon/disp-*` 批量位移 → `03_post` 汇总/绘图；DFPT 仅 2 段）
  ```bash
  vasp phonon -i POSCAR -w ./phonon --supercell 2 2 2 --method disp -j bash          # 仅声子（跳过结构优化）
  vasp relax phonon -i POSCAR -w ./phonon --supercell 2 2 2 --method disp -j bash    # 先做结构优化再做声子
  vasp phonon -i POSCAR -w ./phonon --supercell 2 2 2 --method disp -j bash --submit # 提交当前待办（默认只提交流程首个待办步骤）
  ```
- 分子动力学（1 步：`01_md`，生成 NVT 输入与脚本）
  ```bash
  vasp md -i POSCAR -w ./md --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf          # 准备
  vasp md -i POSCAR -w ./md --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf --submit # 提交
  vasp relax md -i POSCAR -w ./md --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf --submit # 先结构优化再MD
  ```
- 批量并行：目录输入会自动分配子目录，每个结构沿上述流程独立推进
  ```bash
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
