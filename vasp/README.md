# VASP 计算工具 - 使用与开发指引

## 核心概念
- CLI：`vasp/cli.py` 基于 argparse，支持 `--json` 配置合并（优先级：命令行 > JSON > 默认）。`electronic`/`phonon` 已落地，`relax`/`md` 仍为占位。
- Pipeline：`pipelines/base.py` 提供步骤管理、重试、断点续跑（`pipeline_checkpoint.json`）、收敛检测（OUTCAR 关键字），`_submit_job` 仅实现 bash 后台，可扩展 slurm/pbs。
- 具体流程：`pipelines/electronic_properties.py` 自动生成 `01_relax → 02_scf → 03_dos → 04_band → 05_elf → 06_cohp → plots` 并调用 `analysis/plotters.py` 绘图；`pipelines/phonon_properties.py` 负责声子（有限位移/DFPT）。
- 批量：`pipelines/batch.py` 支持目录批量与并行（`--parallel --max-workers`），`pipelines/utils.py` 负责 POTCAR 生成与常用工具。

## 快速上手
- 安装：`pip install -e .`（建议在可写的虚拟环境执行）。
- 查看帮助：`vasp --help`，或 `vasp electronic --help` / `vasp phonon --help`。
- 最小示例：
```bash
vasp electronic -i POSCAR -w ./electronic \
  --kspacing 0.2 --encut 500 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54 --potcar-type PBE

vasp phonon -i POSCAR -w ./phonon --supercell 2 2 2 \
  --method disp --kspacing 0.3 \
  --potcar-dir ~/pot/vasp_pot/potpaw_PBE54 --potcar-type PBE
```
- JSON 示例：`config/config_example.json`，加载方式 `--json config/config_example.json`。

## 任务与目录习惯
- 工作目录会按步骤生成子文件夹（电子流程：`01_relax/02_scf/03_dos/04_band/05_elf/06_cohp/plots`；声子流程类似）。
- 断点文件 `pipeline_checkpoint.json` 支持重复运行续算，`pipeline_report.txt` 汇总各步骤状态。
- 必需自备 POTCAR（`--potcar-dir`、`--potcar-type`），未提供时需手动放置到各子目录。

## 批量计算
- 输入为目录时自动判定批量，或显式 `--batch`。并行控制：`--parallel --max-workers 4`。
- BatchPipeline 会在工作根目录下为每个结构创建子目录，并输出成功/失败计数。

## 开发与扩展
- 新步骤/流程：继承 `BasePipeline`，实现 `get_steps` 与 `execute_step`，在 `cli.py` 增加子命令映射；步骤间数据通过 `steps_data` 传递。
- 新队列系统：在 `_submit_job` 增加 slurm/pbs 提交流程和 job_id 解析；如需自定义完成判据，可覆盖 `_check_job_completed`。
- 新绘图或产物：复用 `analysis/plotters.py`，或在 `pipelines` 中添加新的 `_write_*` 辅助方法。

## 当前限制
- `relax`、`md` 子命令尚未串联真实流程；队列默认仅支持 bash 后台运行。
- 能带高对称路径使用示例值，复杂晶体需要自行调整 KPOINTS。
