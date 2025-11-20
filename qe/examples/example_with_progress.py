#!/usr/bin/env python3
"""
带进度监控的Pipeline执行示例

展示如何在批量计算中使用实时进度监控。

作者：Claude
创建时间：2025-11-20
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.append(str(Path(__file__).parent.parent))

from pipelines.templates.relax_phono_sc import RelaxPhonoSuperconductivityPipeline
from utils.progress import ProgressMonitor
from scheduler.executor import MixedParallelExecutor


def main():
    """主函数：运行带进度监控的批量计算"""

    # 1. 准备结构文件列表
    structure_dir = Path("/home/mayuan/code/my_script/test/H3S")
    structures = list(structure_dir.glob("*.vasp"))

    if not structures:
        print(f"错误: 在 {structure_dir} 中未找到VASP结构文件")
        return

    print(f"找到 {len(structures)} 个结构文件")

    # 2. 配置参数
    config = {
        'ecutwfc': 80,
        'ecutrho': 640,
        'kpoints': (8, 8, 8),
        'qpoints': '6 6 6',
        'pseudopotentials': {
            'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
            'S': 'S.pbe-n-rrkjus_psl.1.0.0.UPF',
        },
        'nprocs': 8,
        'qe_bin': '/home/mayuan/soft/qe-7.4.1/bin/',
    }

    work_dir = Path("./calculations")

    # 3. 创建Pipeline
    pipeline = RelaxPhonoSuperconductivityPipeline(
        structures=structures,
        work_dir=work_dir,
        config=config,
    )

    # 4. 构建任务
    print("\n构建任务列表...")
    pipeline.build_tasks()
    total_tasks = len(pipeline.tasks)

    print(f"总任务数: {total_tasks}")
    print(f"结构数: {len(structures)}")
    print(f"每个结构的步骤数: {total_tasks // len(structures)}")

    # 5. 创建进度监控器
    monitor = ProgressMonitor(total_tasks=total_tasks, refresh_rate=0.5)
    monitor.add_tasks(pipeline.tasks)

    # 6. 创建执行器（带进度回调）
    def progress_callback(task):
        """进度回调函数"""
        monitor.update_task_status(task)

    executor = MixedParallelExecutor(
        max_workers=4,
        progress_callback=progress_callback,
    )

    # 7. 运行计算（带实时监控）
    print("\n开始计算...\n")

    try:
        monitor.run_with_monitor(executor.execute, pipeline.tasks)
    except KeyboardInterrupt:
        print("\n\n计算被用户中断")
        return
    except Exception as e:
        print(f"\n\n计算出错: {e}")
        import traceback
        traceback.print_exc()
        return

    # 8. 打印最终摘要
    monitor.print_summary()

    # 9. 生成报告
    print("\n生成计算报告...")
    pipeline.generate_reports()

    print("\n所有任务完成！")


if __name__ == "__main__":
    main()
