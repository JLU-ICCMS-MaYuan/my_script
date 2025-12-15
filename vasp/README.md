# VASP 计算工具使用指南

## 配置与优先级
- 程序与赝势优先级：① 环境变量 `VASP_STD`/`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`，脚本头可用 `JOB_HEADER_BASH/SLURM/PBS/LSF` 覆盖；② `vasp/config/job_templates.toml` 的 `defaults` 与 `templates`（按 bash/slurm/pbs/lsf 定义头）；③ 兼容旧 `~/.my_scriptrc.py` (`vaspstd_path/vaspgam_path/potcar_dir/*title/default_mpi_procs`)；④ 内置默认。
- 命令行覆盖：`--potcar-dir/--potcar-type`、`-j/--job-system`、`--mpi-procs` 优先于上述配置。`vasp/config_example.json` 仅示例，不会自动读取，只有传 `--json config_example.json` 才会加载。
- 提交方式：`-j/--job-system` 选 `bash/slurm/pbs/lsf`；`--mpi-procs` 可为数字或完整前缀（如 `mpirun -np 16`/`srun -n 16`），未指定取配置或默认 8。

## 目录规则与默认行为
- 结构输入：`-i` 接受单文件或目录，自动判定批量；`--structure-ext` 过滤扩展名（如 `vasp,cif,res,xsf`，默认 vasp）。
- 压强分层：`-p/--pressure` 可多值（单位 GPa，写入 INCAR 的 `PSTRESS=pressure*10`），目录为 `<结构前缀>/<pressure>_GPa/01_relax...`，结构前缀为输入文件去后缀或目录下文件名去后缀。
- 并行批量：目录输入自动批量；`--tasks N` 控制同时处理多少个结构（>1 开启并行），每个压强/结构生成 `batch_summary.txt`。
- 提交策略：默认只生成输入和 `run_*.sh`；加 `--submit` 时提交并等待本命令覆盖的全部步骤完成（批量/多压强全部等待），结束后在工作目录写 `finished`。

## 子命令与步骤依赖
- `vasp relax`：`01_relax`，生成 POSCAR/INCAR/KPOINTS/POTCAR 与 `run_relax.sh`。
- `vasp scf`：`01_relax → 02_scf`，生成 CHGCAR。
- `vasp dos|band|elf|cohp|bader|fermisurface`：默认 `01_relax → 02_scf → <目标步骤>`（Bader/费米面目录分别为 `07_bader`、`08_fermi`）。
- `vasp phonon`：固定包含结构优化，流程 `01_relax → 02_phonon (disp/dfpt) → phonon_band → phonon_dos → plots`。
- `vasp md`：固定 `01_relax → 02_md`，自动继承 CONTCAR/POTCAR。
- 所有命令在 `bash` 下直接 `mpirun -np <N> vasp_std`；在队列模式下 `sbatch/qsub/bsub run_*.sh`，任务号消失后再检查 OUTCAR 判定完成。

## 常用示例（均会自动建目录）
- 单结构：
  ```bash
  vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs 48 --submit
  vasp dos -i POSCAR --submit
  vasp phonon -i POSCAR --supercell 2 2 2 --method disp --submit
  vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 --submit
  ```
- 目录批量 + 并行：
  ```bash
  vasp relax -i ./structures --structure-ext vasp,cif --pressure 0 10 --tasks 4 -j slurm
  vasp band  -i ./structures --tasks 3 --submit
  ```
- 多模块组合（自动补齐依赖，relax 完成后并行调度 phonon 与电子流程）：
  ```bash
  vasp combo relax phonon dos -i ./stdlibs/ -p 0 5 -j slurm --encut 600 --kspacing 0.18 --mpi-procs "srun -n 32" --submit
  vasp combo relax md -i POSCAR --potim 1.0 --nsw 200 --submit
  ```

## 输出与断点
- 每步写 `pipeline_checkpoint.json` 支持续跑，`pipeline_report.txt` 记录状态；成功结束写 `finished`。
- 批量模式在压强目录下生成 `batch_summary.txt` 汇总成功/失败。***


## 测试所有命令行功能

下面给出覆盖全部功能的命令示例（均会按 <结构前缀>/<pressure>_GPa/01_relax... 建目录），默认不加 --submit 仅准备输入；加 --submit 会提
  交并等待完成。

  - 基础结构优化
    vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs "mpirun -np 32" --submit
  - 自洽（带 relax）
    vasp scf -i POSCAR --kspacing 0.2 --encut 520 -j bash --mpi-procs 16 --submit
  - DOS / Band / ELF / COHP / Bader / Fermi surface（都自动补齐 relax→scf）
    vasp dos  -i POSCAR --encut 550 --kspacing 0.18 --mpi-procs "srun -n 48" --submit
    vasp band -i POSCAR --mpi-procs 32 --submit
    vasp elf  -i POSCAR --submit
    vasp cohp -i POSCAR --submit
    vasp bader -i POSCAR --submit
    vasp fermisurface -i POSCAR --submit
  - 声子（disp/dfpt）
    vasp phonon -i POSCAR --supercell 2 2 2 --method disp -j bash --mpi-procs "mpirun -np 24" --submit
    vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 500 -j slurm --mpi-procs 32 --submit
  - 多模块组合（relax 完成后并行跑 phonon 与电子链路）
    vasp combo relax phonon dos -i POSCAR -p 0 5 -j slurm --encut 600 --kspacing 0.18 --mpi-procs "srun -n 32" --submit
  - 批量 + 并行结构（目录输入，扩展过滤）
    vasp relax -i ./structures --structure-ext cif,xsf -p 0 10 --tasks 4 -j slurm --submit
    vasp band  -i ./structures --tasks 3 -j bash --mpi-procs "mpirun -np 8" --submit
  - 多压强并行示例（目录批量）
    vasp combo relax scf dos -i ./stdlibs/ -p 0 50 100 --tasks 3 -j slurm --mpi-procs "mpirun -np 16" --submit
  - 自定义启动器字符串示例
    vasp scf -i POSCAR -j slurm --mpi-procs "ibrun -n 64" --submit
  - 仅准备、不提交示例（检查生成的 run_*.sh）
    vasp combo relax phonon dos -i POSCAR --kspacing 0.2 --encut 520 -j bash --mpi-procs 8

  这些命令覆盖：单结构/批量、多压强、所有模块、不同队列/启动器字符串、并行 tasks、submit/prepare-only。
