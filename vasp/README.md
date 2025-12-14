# VASP 计算工具 - 使用与开发指引

## 配置来源
- 程序路径与赝势：优先读取环境变量 `VASP_STD`、`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`；其次读取 `vasp/config/job_templates.toml` 的 defaults；最后回退到内置默认（兼容旧 `.my_scriptrc.py`）。
- 队列脚本头：`vasp/config/job_templates.toml` 中按 bash/slurm/pbs/lsf 提供模板，可用环境变量 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖。
- 队列选择：`-j/--job-system` 取 `bash/slurm/pbs/lsf`；MPI 进程数可用 `--mpi-procs` 指定（默认 8 或配置值）。
- POTCAR：默认使用 `POTCAR_DIR` 或模板 defaults，可被 `--potcar-dir` 覆盖；`--potcar-type` 选择 PBE/LDA/PW91。

## 执行模式（统一约定）
- 默认仅准备：不加 `--submit` 时只生成当前待执行步骤的输入与 `run_*.sh`，不提交、不等待；多步流程需多次调用推进下一步（断点文件自动跳过已完成步骤）。
- 提交并等待：加 `--submit` 时提交并等待当前命令覆盖的全部步骤完成（含子作业/批量/多压强所有结构）；完成后工作目录写入 `finished` 标记。
- Bash/队列：`-j bash` 直接在目录下 `mpirun -np <N> vasp_std` 前台跑；`-j slurm/pbs/lsf` 则 `sbatch/qsub/bsub` 脚本并轮询队列，任务号消失后再检查 OUTCAR 判定完成；`--mpi-procs` 控制 `-np`。
- 外压与目录：`-p/--pressure` 支持单/多压强（单位 GPa，写入 INCAR 的 `PSTRESS`，kBar），工作目录自动分层为 `<pressure>_GPa/01_relax...`。
- 结构输入：`-i` 可为单文件或目录，自动判定批量；`--structure-ext` 过滤目录下的结构类型（如 `vasp`/`cif`/`res`，默认 vasp）。
  - 目录输入示例：`vasp relax -i stdlibs/ -j slurm --mpi-procs 8` 会扫描 `stdlibs/` 下匹配扩展（默认 vasp，可用 `--structure-ext cif,res,xsf`），并生成 `stdlibs/0_GPa/<结构名>/01_relax...`。

## 各命令的执行流程
- 结构优化 `vasp relax`：生成 `01_relax`（POSCAR/INCAR/KPOINTS/POTCAR + run_relax.sh）→ 提交/准备。工作根目录自动取输入文件名去后缀（目录输入则用目录名）并在其下按压强分层。
- 支持目录输入与多压强：`-i ./structures --structure-ext cif -p 0 5` 会按压强/结构名分层生成计算目录。
- SCF `vasp scf`：默认步骤 `01_relax → 02_scf`。
- DOS/Band/ELF/COHP/Bader/Fermi surface：分别使用 `vasp dos|band|elf|cohp|bader|fermi`，默认都会执行 `01_relax → 02_scf → <目标步骤>`（bader/fermi 使用 scf 结果）；可通过 `--steps` 调整子集，但需保证 relax 和 scf 在前。
- 声子 `vasp phonon`：
  - 默认仅声子：`02_phonon` 下生成位移超胞（disp）或 DFPT 输入 → 对每个 `disp-*` 准备/提交 → phonopy 后处理能带/DOS → 绘图。
  - 若需先优化：用 `vasp relax phonon ...` 或 `vasp phonon --with-relax ...`，流程变为 `01_relax` 后再进入声子步骤。
- 分子动力学 `vasp md`：
  - 默认仅 MD：生成 `01_md` 输入与脚本。
  - 若需先优化：`vasp relax md ...` 或 `vasp md --with-relax ...`，先 `01_relax`，再 `02_md`，并自动继承 CONTCAR/POTCAR。
- 批量模式：输入为目录或显式 `--batch` 时，扫描结构文件，逐个创建子目录独立运行；`--tasks N` 控制同时运行的结构数（默认串行）。
  - 目录层级：`<cwd>/<结构名>/<pressure>_GPa/01_relax...`；单文件同样使用 `<cwd>/<文件名去后缀>/<pressure>_GPa/...`。

## 快速示例与步骤说明
- 结构优化（1 步：`01_relax`，准备 INCAR/KPOINTS/POTCAR + `run_relax.sh`）
  ```bash
  vasp relax -i POSCAR -j slurm --mpi-procs 48          # 仅准备，工作目录为 ./POSCAR/0_GPa/01_relax
  vasp relax -i POSCAR -j slurm --mpi-procs 48 --submit # 提交并等待
  vasp relax -i POSCAR -p 0 5 10                        # 0/5/10 GPa 三组输入，各自生成 <P>_GPa 子目录
  ```
- 电子性质（6~7 步：`01_relax → 02_scf → 03_dos → 04_band → [05_elf] → [06_cohp] → plots`）
  ```bash
  vasp electronic -i POSCAR --kspacing 0.2 --encut 500 \
    --include-elf --include-cohp --potcar-type PBE -j pbs                # 准备当前步
  vasp electronic -i POSCAR --kspacing 0.2 --encut 500 \
    --include-elf --include-cohp --potcar-type PBE -j pbs --submit       # 提交并等待
  vasp electronic -i POSCAR -p 0 5 --kspacing 0.2                        # 0/5 GPa 各生成子目录（自动以文件名为根目录）
  # 提示：每次运行只处理一个待办步骤；当 relax 完成后再运行同一命令推进到 scf，依次推进到 band/elf/cohp，最后自动绘图。
  ```
- 声子（位移法典型 3 段：`01_relax → 02_phonon/disp-*` 批量位移 → `03_post` 汇总/绘图；DFPT 仅 2 段）
  ```bash
  vasp phonon -i POSCAR --supercell 2 2 2 --method disp -j bash          # 仅声子（跳过结构优化）
  vasp relax phonon -i POSCAR --supercell 2 2 2 --method disp -j bash    # 先做结构优化再做声子
  vasp phonon -i POSCAR --supercell 2 2 2 --method disp -j bash --submit # 提交并等待
  ```
- 分子动力学（1 步：`01_md`，生成 NVT 输入与脚本）
  ```bash
  vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf          # 准备
  vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf --submit # 提交并等待
  vasp relax md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 -j lsf --submit # 先结构优化再MD
  ```
- 批量并行：目录输入会自动分配子目录，每个结构沿上述流程独立推进
  ```bash
  vasp electronic -i ./structures --tasks 4
  vasp relax -i ./structures --structure-ext cif -p 0 10 --tasks 8
  ```

## 批量与并行
- 目录输入或 `--batch` 触发批量模式；`--tasks N` 控制并行。
- 每个结构生成独立子目录并统计成功/失败，汇总写入 `batch_summary.txt`。

## 工作目录与产物
- 电子流程：分步目录与 `plots/`（band/dos/elf/cohp）；
- 声子流程：`02_phonon/disp-*` + phonopy 后处理产物与 `plots/`；
- 断点文件 `pipeline_checkpoint.json` 支持续跑，报告 `pipeline_report.txt` 便于回溯；成功结束会生成 `finished`。

## 注意事项
- 队列提交通过 rc/默认脚本头生成并调用 `mpirun`；命令缺失时回退本地 bash 执行，bash 环境即可直接跑。
- `_wait_for_job` 同时检查队列状态（slurm/pbs/lsf）与 OUTCAR 关键字；队列号消失会再检查完成标志后才报错。
- 仓库不内置 VASP 可执行文件/赝势，请确保路径、权限与 MPI 环境正确。
