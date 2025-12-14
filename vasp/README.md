# VASP 计算工具使用说明

## 配置优先级
- 程序与赝势：优先读环境变量 `VASP_STD`/`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`；其次读 `vasp/config/job_templates.toml` defaults；最后回退旧 `.my_scriptrc.py`。
- 队列脚本头：`vasp/config/job_templates.toml` 中按 bash/slurm/pbs/lsf 定义，可用环境变量 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖。
- 任务提交：`-j/--job-system` 选 bash/slurm/pbs/lsf，`--mpi-procs` 控制 mpirun -np，未指定时取配置或默认 8。

## 目录与默认行为
- 结构输入：`-i` 接受单文件或目录，自动判定批量；`--structure-ext` 过滤扩展（如 `vasp,cif,res,xsf`），默认 vasp。
- 压强层级：`-p/--pressure` 支持多值（GPa），写入 INCAR 的 `PSTRESS`（kBar）。工作目录统一为 `<cwd>/<结构前缀>/<pressure>_GPa/01_relax...`。
- 默认只准备：不加 `--submit` 时仅生成输入与 `run_*.sh`，写断点文件方便续跑；加 `--submit` 时提交并轮询队列/OUTCAR，流程全部结束后写 `finished`。
- 并行与批量：目录输入自动批量；`--tasks N` 控制同时处理的结构数（N>1 并行）。每个结构/压强单独生成 `batch_summary.txt` 汇总。

## 核心命令与步骤
- `vasp relax`：1 步 → `01_relax`（POSCAR/INCAR/KPOINTS/POTCAR + run_relax.sh）。  
  示例：`vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs 48 --submit`
- `vasp scf`：2 步 → `01_relax → 02_scf`。scf 使用 relax 结构并生成 CHGCAR。  
  示例：`vasp scf -i POSCAR --kspacing 0.2 --submit`
- `vasp dos|band|elf|cohp|bader|fermisurface`：默认 `01_relax → 02_scf → <目标步骤>`（bader/fermi 目录为 `07_bader`、`08_fermi`）。  
  示例：`vasp dos -i POSCAR --kspacing 0.2 --steps relax,scf,dos,band --submit`（可通过 `--steps` 拼装如 relax→scf→dos→band→elf）。
- `vasp phonon`：默认包含结构优化，流程 `01_relax → 02_phonon (disp/dfpt) → phonon_band/phonon_dos → plots`。  
  示例：`vasp phonon -i POSCAR --supercell 2 2 2 --method disp --submit`
- `vasp md`：固定 `01_relax → 02_md`，自动继承 CONTCAR/POTCAR。  
  示例：`vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 --submit`
- 依赖约束：dos/band/elf/cohp/bader/fermi 都依赖 relax+scf；phonon/MD 必须在 relax 之后（不依赖 scf）；`--steps` 会自动补齐依赖顺序。

## 批量与多压强示例
- 多结构 + 多压强 + 并行：  
  `vasp relax -i ./stdlibs --structure-ext vasp,cif --pressure 0 5 10 --tasks 4 -j slurm`  
  生成 `<cwd>/<结构名>/<P>_GPa/01_relax...`，完成后各目录写 `finished`，汇总于 `<结构名>/<P>_GPa/batch_summary.txt`。
- 串联示例（积木式组合）：  
  `vasp dos -i POSCAR --steps relax,scf,dos,band,elf,cohp,bader,fermisurface --submit` 一次跑完全部电子性质链路；  
  `vasp phonon -i POSCAR --supercell 3 3 3 --method dfpt --submit` 跑完整声子流程。

## 常用参数速查
- `-i/--input` 结构文件或目录；`--structure-ext` 过滤扩展名。
- `-p/--pressure` 多压强（GPa）；目录命名 `<pressure>_GPa`，INCAR 写 `PSTRESS=pressure*10`。
- `-j/--job-system` bash/slurm/pbs/lsf；`--mpi-procs` 进程数。
- `--tasks` 并行结构数（>1 启动并行）；`--steps` 自定义步骤序列（自动插入 relax/scf 依赖）。
- `--supercell`、`--method`（phonon）；`--potim --tebeg --teend --nsw`（md）；`--kspacing --encut --dos-type`（电子性质相关）。
