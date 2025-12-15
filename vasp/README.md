# VASP 计算工具使用指南

## 配置与优先级
- 程序与赝势来源（高→低）：① 环境变量 `VASP_STD`/`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`，脚本头可用 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖；② `vasp/config/job_templates.toml` 的 `defaults` 与 `templates`（按 bash/slurm/pbs/lsf 定义头）；③ 兼容旧 `~/.my_scriptrc.py` (`vaspstd_path/vaspgam_path/potcar_dir/*title/default_mpi_procs`)；④ 内置默认。
- 命令行优先：`--potcar-dir/--potcar-type`、`-j/--job-system`、`--mpi-procs`、`--encut`、`--kspacing` 等最高优先级。`vasp/config_example.json` 仅示例，不会自动读取，只有传 `--json config_example.json` 才会加载。
- 提交方式：`-j/--job-system` 选 `bash/slurm/pbs/lsf`；`--mpi-procs` 可为数字或完整前缀（如 `mpirun -np 16`、`srun -n 16`），未指定取配置或默认 8。
- 赝势管理：优先使用当前工作目录下的 `potcar_lib`（支持 `potcar_lib/元素` 或 `potcar_lib/元素/POTCAR`），若缺元素且提供 `potcar_dir` 则仅在 `potcar_dir`（含 `potcar_type` 子目录）中寻找唯一候选并复制到 `potcar_lib`。若找不到或候选不唯一会报错，请手动将所需 POTCAR 放入 `potcar_lib` 后重试。

## 目录规则与执行模式
- 结构输入：`-i` 接受单文件或目录，自动判定批量；`--structure-ext` 过滤扩展名（如 `vasp,cif,res,xsf`，默认 vasp）。
- 压强分层：`-p/--pressure` 支持多值（GPa，写入 INCAR `PSTRESS=pressure*10`），目录为 `<结构前缀>/<pressure>_GPa/01_relax...`，结构前缀为文件去后缀或目录下文件名去后缀。
- 批量与并行：目录输入自动批量；`--tasks N` 控制同时处理的结构数（>1 并行），每个压强/结构生成 `batch_summary.txt`。
- 提交策略：默认仅生成输入与 `run_*.sh`；加 `--submit` 时提交并等待本命令覆盖的全部步骤完成（含批量/多压强、phonon 子位移），结束后在工作目录写 `finished`。
- 状态与重跑：`prepare_only` 模式步骤标记为 `PREPARED`，不会记为完成；重复运行会按产物检查（OUTCAR/CHGCAR 等）决定是否重跑，缺产物或仅 `PREPARED` 一律视为未完成。失败不重试，直接报错退出。

## 步骤依赖与并行编排
- 基础步骤：`relax`；`scf` 依赖 relax；`dos/band/elf/cohp/bader/fermisurface` 依赖 relax+scf；`phonon` 依赖 relax；`md` 依赖 relax。
- 组合执行：`combo` 支持串写模块（如 `relax phonon dos`），程序自动补齐依赖；relax 完成后并行启动 phonon 与电子链路（scf→dos/band/elf/cohp/bader/fermisurface），md 也可并行，统一由 Python 监控。
- 任务提交：bash 直接在目录内启动；队列模式通过 `run_*.sh` 调用 `sbatch/qsub/bsub`，任务号消失后再查 OUTCAR 判定完成。

## 输出、断点与检查
- 每条子链有独立 `pipeline_checkpoint.json` 与 `pipeline_report.txt`；成功结束写 `finished`。已完成且产物存在才会跳过，若缺少关键文件则降级为未完成重跑。
- 批量模式在压强目录下生成 `batch_summary.txt` 汇总成功/失败。

## 常用示例（自动建目录）
- 单结构：
  ```bash
  vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs 48 --submit
  vasp dos   -i POSCAR --submit
  vasp phonon -i POSCAR --supercell 2 2 2 --method disp --submit
  vasp md    -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 --submit
  ```
- 目录批量 + 并行：
  ```bash
  vasp relax -i ./structures --structure-ext vasp,cif --pressure 0 10 --tasks 4 -j slurm
  vasp band  -i ./structures --tasks 3 --submit
  ```
- 多模块组合（自动补齐依赖，relax 后并行调度 phonon 与电子链路）：
  ```bash
  vasp combo relax phonon dos -i ./stdlibs/ -p 0 5 -j slurm --encut 600 --kspacing 0.18 --mpi-procs "srun -n 32" --submit
  vasp combo relax md -i POSCAR --potim 1.0 --nsw 200 --submit
  ```
- 运行前准备赝势：在当前工作目录创建 `potcar_lib`，放入选定赝势（如 `potcar_lib/Si/POTCAR`）；缺失元素时若提供 `potcar_dir` 会复制唯一匹配项到 `potcar_lib`，否则报错提示手动补齐。
- 全功能回归示例（覆盖全部子命令与参数形态）：
  - 基础优化：`vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs "mpirun -np 32" --submit`
  - 自洽：`vasp scf -i POSCAR --kspacing 0.2 --encut 520 -j bash --mpi-procs 16 --submit`
  - DOS/Band/ELF/COHP/Bader/Fermi：`vasp dos -i POSCAR --encut 550 --kspacing 0.18 --mpi-procs "srun -n 48" --submit`；`vasp band -i POSCAR --mpi-procs 32 --submit`；`vasp elf|cohp|bader|fermisurface -i POSCAR --submit`
  - 声子：`vasp phonon -i POSCAR --supercell 2 2 2 --method disp -j bash --mpi-procs "mpirun -np 24" --submit`
  - MD：`vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 500 -j slurm --mpi-procs 32 --submit`
  - 批量+多压强：`vasp combo relax scf dos -i ./stdlibs/ -p 0 50 100 --tasks 3 -j slurm --mpi-procs "mpirun -np 16" --submit`
  - 自定义启动器：`vasp scf -i POSCAR -j slurm --mpi-procs "ibrun -n 64" --submit`
  - 仅准备：`vasp combo relax phonon dos -i POSCAR --kspacing 0.2 --encut 520 -j bash --mpi-procs 8`

这些命令覆盖单结构/批量、多压强、全部模块、不同队列/启动器字符串、并行 tasks，以及 submit 与 prepare-only 场景。***
