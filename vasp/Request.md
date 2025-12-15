# VASP 需求总览（持续维护）

## 运行行为与依赖
- 子命令：relax、scf、dos、band、elf、cohp、bader、fermisurface、phonon、md，支持 `combo` 串写（如 `relax phonon dos`），自动补齐依赖。电子链路依赖 relax+scf，phonon/md 依赖 relax，phonon 与 scf 链路可在 relax 后并行。
- 默认仅准备输入与脚本；`--submit` 时提交并等待本命令覆盖的全部步骤（批量/多压强/phonon 子位移）。
- prepare_only 步骤标记为 `PREPARED`，不算完成；缺产物或 PREPARED 一律视为未完成可重跑。失败不重试，直接报错并停在当前步骤。
- 目录/压强规则：工作目录 `<结构前缀>/<pressure>_GPa/01_relax...`，结构前缀取文件去后缀或目录下文件名去后缀。批量自动扫描目录；`--tasks` 控制并行结构数。

## 队列与提交
- `--mpi-procs` 可为数字或完整启动前缀（如 `mpirun -np 16`、`srun -n 16`），按脚本所在目录执行 sbatch/qsub/bsub；提交失败直接报错，不回退本地。
- Slurm/LSF 解析输出提取纯数字 job_id；队列轮询结合 OUTCAR 完成判定。

## POTCAR 管理
- 不再自动选择赝势。优先使用当前工作目录下 `potcar_lib`（支持 `potcar_lib/元素` 或 `potcar_lib/元素/POTCAR`），按 POSCAR 第 6 行元素顺序合并。
- 若缺元素且提供 `potcar_dir`，仅在 `potcar_dir`（含 `potcar_type` 子目录）中寻找唯一候选并复制到 `potcar_lib`；找不到或多候选则报错，需手动放入 `potcar_lib` 后重试。

## 断点与文件
- checkpoint 仅保留 `steps_status` 与 `steps_data`（传递 relaxed_structure、CHGCAR 等路径）；无用的 timestamp 已移除。
- 每子链生成独立报告与 `finished` 标记；批量汇总 `batch_summary.txt`。

## 已知约束/待关注
- phonon/band 等后处理在 prepare-only 下不应假定产出文件（如 band.yaml、vasprun.xml）；已按失败直接停的策略处理。
- 需要保持 README 与行为同步更新，新增需求按模块分类增补本文件而非简单追加。***
