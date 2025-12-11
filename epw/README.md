# EPW 电声耦合工具 - 使用与开发指引

## 状态概览
- CLI 位于 `epw/cli.py`，子命令统一为 `epw run`，支持 `--mode full/eband/phono/elph/sc` 与 `--json` 配置合并（命令行优先）。当前逻辑仅打印配置，具体 Workflow 尚未接入。
- 依赖 QE 的 SCF/Phonon 输出与 Wannier90，可用于电声耦合、超导与输运研究。

## 快速上手（配置预检）
- 安装：`pip install -e .`。
- 常用命令（当前只校验参数，不提交计算）：
```bash
epw run -i scf.out -w ./epw_full --mode full \
  --nkf 20 20 20 --nqf 20 20 20 --pp-dir ~/pot/qe_pp

epw run -i scf.out -w ./epw_elph --mode elph \
  --nk 8 8 8 --nq 4 4 4 --nkf 16 16 16 --nqf 8 8 8
```
- JSON 示例：`config_example.json`，通过 `--json config_example.json` 载入。

## 代码结构速览
- `cli.py`：参数解析与配置合并，负责根据 `mode` 打印当前任务设定。
- `analysis/`、`utils/`：预留分析与工具函数目录，可复用 QE 输出解析、Wannier 后处理。
- `config/` 与示例 JSON：提供默认网格、投影、展宽等字段，便于批量生成配置。

## 开发与扩展指引
- Workflow 设计：建议新增 `epw/workflows.py`（当前在 `cli.py` 中被 TODO 注释），按 `mode` 划分执行链，如 NSCF → Wannier 拟合 → EPW → 超导/输运分析，并写入中间文件与报告。
- 队列与监控：可参考 VASP Pipeline 的 `_submit_job`/`_check_job_completed` 思路，实现 bash/slurm/pbs 提交与轮询。
- 数据依赖：在 Workflow 中显式校验 QE/Phonon/Wannier 文件路径（如 DVSCF、amn/mmn/eig），缺失时提前失败并提示。
- 结果复用：将关键路径写入工作目录下的 meta/JSON，方便后续 elph/sc 模式继续计算。

## 使用提示
- 确保 `epw.x`、`pw.x`、`ph.x`、`wannier90.x` 在 PATH 中；`--pp-dir` 仅传递赝势位置，仓库不会自动下载。
- 电声耦合计算对网格与展宽敏感，建议在 JSON 中固定 nk/nq/nkf/nqf 与 `fsthick/degaussw/degaussq`，避免命令行遗漏。
