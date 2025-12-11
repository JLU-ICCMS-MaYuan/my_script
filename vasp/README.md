# VASP 计算工具 - 使用与开发指引

## 功能概览
- CLI：`vasp/cli.py` 提供 `electronic`、`phonon`、`relax`、`md` 子命令，支持 `--json` 配置合并（优先级：命令行 > JSON > 默认）。
- Pipeline：`pipelines/base.py` 负责步骤管理、重试、断点续跑（`pipeline_checkpoint.json`）、收敛检测（OUTCAR 关键字）；提交脚本默认仅支持 bash 后台，可扩展 slurm/pbs。
- 电子流程：`pipelines/electronic_properties.py` 自动生成 `01_relax → 02_scf → 03_dos → 04_band → 05_elf → 06_cohp → plots`，并调用 `analysis/plotters.py` 绘图。
- 声子流程：`pipelines/phonon_properties.py` 支持有限位移(`disp`)与 DFPT(`dfpt`)；批量与并行由 `pipelines/batch.py` 提供。
- 状态：`electronic/phonon` 已落地；`relax/md` 仍为占位，需要补充 Workflow 后才能实际提交。

## 前置依赖
- VASP 可执行文件（默认脚本使用 `mpirun -np 8 ~/soft/vasp.6.3.2/bin/vasp_std`，需按环境修改），具备 MPI 运行权限。
- POTCAR 资源：通过 `--potcar-dir`、`--potcar-type` 指定；未提供时需手动放置到各步骤目录。
- 推荐在可写虚拟环境中 `pip install -e .`，并准备 `config/config_example.json` 作为模板。

## 快速使用
```bash
# 电子性质全流程（含 ELF/COHP，单结构）
vasp electronic -i POSCAR -w ./electronic \
  --kspacing 0.2 --encut 500 \
  --include-elf --include-cohp \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54 --potcar-type PBE

# 声子：有限位移超胞 2×2×2
vasp phonon -i POSCAR -w ./phonon \
  --supercell 2 2 2 --method disp --kspacing 0.3 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54 --potcar-type PBE

# 批量/并行（目录自动识别结构文件）
vasp electronic -i ./structures -w ./batch_elec --parallel --max-workers 4 \
  --kspacing 0.25 --potcar-dir ~/pot/vasp_pot/potpaw_PBE54

# JSON 配置
vasp phonon -i POSCAR -w ./phonon --json config/config_example.json
```

## 批量与并行
- 输入为目录时自动判定批量，也可显式 `--batch`；并行控制：`--parallel --max-workers N`。
- BatchPipeline 会在工作根目录下为每个结构创建独立子目录，并统计成功/失败。

## 工作目录与产物
- 电子流程生成分步子目录及 `plots/` 图像；声子流程生成对应步骤目录。
- 断点：`pipeline_checkpoint.json` 支持重复运行续算；`pipeline_report.txt` 汇总步骤状态。
- KPOINTS 默认采用自动/示例路径，如需精确高对称线需自定义。

## 配置与优先级
- 配置合并顺序：命令行 > `--json` > 代码默认值；建议将通用参数写入 JSON，差异参数用命令行覆盖。
- 可通过 `--log-level` 控制日志；未提供 `--potcar-dir` 时需提前准备 POTCAR。

## 局限与注意
- `relax`、`md` 未接入真实 Workflow；队列仅实现 bash 提交，slurm/pbs 需自行扩展 `_submit_job`。
- 实际计算依赖外部 VASP 与队列环境，仓库不内置可执行文件，运行前需自行校验路径、MPI 及文件读写权限。
