# VASP 计算工具使用指南

## 配置与优先级
- 程序与赝势：优先读取环境变量 `VASP_STD`/`VASP_GAM`、`POTCAR_DIR`、`VASP_MPI_PROCS`；其次读取 `vasp/config/job_templates.toml` 中的 `defaults`；最后兼容旧 `.my_scriptrc.py`。
- 队列脚本头：`job_templates.toml` 中按 `bash/slurm/pbs/lsf` 提供模板，可用 `JOB_HEADER_BASH/SLURM/PBS/LSF` 环境变量覆盖。
- 提交方式：`-j/--job-system` 选 `bash/slurm/pbs/lsf`；`--mpi-procs` 控制 `mpirun -np`，未指定取配置或默认 8。

## 目录规则与默认行为
- 结构输入：`-i` 接受单文件或目录，自动判定批量；`--structure-ext` 过滤扩展名（如 `vasp,cif,res,xsf`，默认 vasp）。
- 压强分层：`-p/--pressure` 可多值（单位 GPa，写入 INCAR 的 `PSTRESS=pressure*10`），目录为 `<结构前缀>/<pressure>_GPa/01_relax...`，结构前缀为输入文件去后缀或目录下文件名去后缀。
- 并行批量：目录输入自动批量；`--tasks N` 控制同时处理多少个结构（>1 开启并行），每个压强/结构生成 `batch_summary.txt`。
- 提交策略：默认只生成输入和 `run_*.sh`；加 `--submit` 时提交并等待本命令覆盖的全部步骤完成（批量/多压强全部等待），结束后在工作目录写 `finished`。

## 子命令与步骤依赖
- `vasp relax`：`01_relax`，生成 POSCAR/INCAR/KPOINTS/POTCAR 与 `run_relax.sh`。
- `vasp scf`：`01_relax → 02_scf`，生成 CHGCAR。
- `vasp dos|band|elf|cohp|bader|fermisurface`：默认 `01_relax → 02_scf → <目标步骤>`（Bader/费米面目录分别为 `07_bader`、`08_fermi`）。`--steps` 可自定义序列，程序自动补齐依赖（确保 relax、scf 在前）。
- `vasp phonon`：固定包含结构优化，流程 `01_relax → 02_phonon (disp/dfpt) → phonon_band → phonon_dos → plots`。
- `vasp md`：固定 `01_relax → 02_md`，自动继承 CONTCAR/POTCAR。
- 所有命令在 `bash` 下直接 `mpirun -np <N> vasp_std`；在队列模式下 `sbatch/qsub/bsub run_*.sh`，任务号消失后再检查 OUTCAR 判定完成。

## 常用示例（均会自动建目录）
- 单结构：
  ```bash
  vasp relax -i POSCAR -p 0 5 -j slurm --mpi-procs 48 --submit
  vasp dos -i POSCAR --steps relax,scf,dos,band,elf,cohp,bader,fermisurface --submit
  vasp phonon -i POSCAR --supercell 2 2 2 --method disp --submit
  vasp md -i POSCAR --potim 1.0 --tebeg 300 --teend 300 --nsw 200 --submit
  ```
- 目录批量 + 并行：
  ```bash
  vasp relax -i ./structures --structure-ext vasp,cif --pressure 0 10 --tasks 4 -j slurm
  vasp band  -i ./structures --tasks 3 --submit
  ```

## 输出与断点
- 每步写 `pipeline_checkpoint.json` 支持续跑，`pipeline_report.txt` 记录状态；成功结束写 `finished`。
- 批量模式在压强目录下生成 `batch_summary.txt` 汇总成功/失败。***
