"""
VASP Pipeline批量计算

支持批量处理多个结构文件

作者：Claude
创建时间：2025-11-20
"""

import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Type
from concurrent.futures import ProcessPoolExecutor, as_completed

from vasp.pipelines.base import BasePipeline
from vasp.pipelines.utils import validate_structure_file, generate_summary_report

logger = logging.getLogger(__name__)


class BatchPipeline:
    """
    批量Pipeline执行器

    读取目录中的所有.vasp文件，批量执行指定的Pipeline
    支持串行和并行两种模式
    """

    def __init__(
        self,
        pipeline_class: Type[BasePipeline],
        structures_dir: Path,
        work_root: Path,
        pipeline_kwargs: Optional[Dict[str, Any]] = None,
        parallel: bool = False,
        max_workers: int = 4,
        structure_exts: Optional[List[str]] = None,
        pressure_label: Optional[str] = None,
    ):
        """
        初始化批量Pipeline

        Parameters
        ----------
        pipeline_class : Type[BasePipeline]
            要执行的Pipeline类（如PropertiesPipeline）
        structures_dir : Path
            包含结构文件的目录
        work_root : Path
            工作根目录，会为每个结构创建子目录
        pipeline_kwargs : Dict, optional
            传递给Pipeline的额外参数
        parallel : bool
            是否并行执行（默认串行）
        max_workers : int
            并行执行时的最大worker数
        """
        self.pipeline_class = pipeline_class
        self.structures_dir = Path(structures_dir)
        self.work_root = Path(work_root)
        self.pipeline_kwargs = pipeline_kwargs or {}
        self.parallel = parallel
        self.max_workers = max_workers
        self.structure_exts = structure_exts or ["vasp"]
        self.pressure_label = pressure_label

        # 创建工作根目录
        self.work_root.mkdir(parents=True, exist_ok=True)

        # 扫描结构文件
        self.structure_files = self._scan_structure_files()

        logger.info(f"找到 {len(self.structure_files)} 个结构文件")

    def _scan_structure_files(self) -> List[Path]:
        """扫描目录中的结构文件"""
        if not self.structures_dir.exists():
            raise FileNotFoundError(f"结构目录不存在: {self.structures_dir}")

        structure_files = []

        patterns = []
        for ext in self.structure_exts:
            e = ext.lower()
            if e == "vasp":
                patterns.extend(['*.vasp', '*.POSCAR', 'POSCAR*'])
            elif e == "cif":
                patterns.append("*.cif")
            elif e == "res":
                patterns.append("*.res")
            elif e == "xsf":
                patterns.append("*.xsf")
            else:
                patterns.append(f"*.{e}")

        for pattern in patterns:
            structure_files.extend(self.structures_dir.glob(pattern))

        # 去重并验证
        valid_files = []
        seen = set()

        for f in structure_files:
            if f in seen:
                continue

            if validate_structure_file(f):
                valid_files.append(f)
                seen.add(f)
            else:
                logger.warning(f"跳过无效结构文件: {f}")

        return sorted(valid_files)

    def run(self) -> List[Dict[str, Any]]:
        """
        执行批量计算

        Returns
        -------
        List[Dict]
            每个结构的计算结果
        """
        logger.info("="*60)
        logger.info(f"开始批量计算: {self.pipeline_class.__name__}")
        logger.info(f"结构数量: {len(self.structure_files)}")
        logger.info(f"执行模式: {'并行' if self.parallel else '串行'}")
        logger.info("="*60 + "\n")

        results = []

        if self.parallel:
            results = self._run_parallel()
        else:
            results = self._run_serial()

        # 生成汇总报告
        summary_file = self.work_root / "batch_summary.txt"
        generate_summary_report(results, summary_file)

        # 统计
        total = len(results)
        success = sum(1 for r in results if r.get('success'))
        failed = total - success

        logger.info("\n" + "="*60)
        logger.info("批量计算完成！")
        logger.info(f"总计: {total}, 成功: {success}, 失败: {failed}")
        logger.info(f"汇总报告: {summary_file}")
        logger.info("="*60)

        return results

    def _run_serial(self) -> List[Dict[str, Any]]:
        """串行执行"""
        results = []

        for i, structure_file in enumerate(self.structure_files, 1):
            logger.info(f"\n处理结构 {i}/{len(self.structure_files)}: {structure_file.name}")

            result = self._run_single_structure(structure_file)
            results.append(result)

            if not result.get('success'):
                logger.warning(f"结构 {structure_file.name} 计算失败，继续下一个")

        return results

    def _run_parallel(self) -> List[Dict[str, Any]]:
        """并行执行"""
        results = []

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # 提交所有任务
            futures = {
                executor.submit(self._run_single_structure, structure_file): structure_file
                for structure_file in self.structure_files
            }

            # 收集结果
            for future in as_completed(futures):
                structure_file = futures[future]

                try:
                    result = future.result()
                    results.append(result)

                    status = "成功" if result.get('success') else "失败"
                    logger.info(f"结构 {structure_file.name}: {status}")

                except Exception as e:
                    logger.error(f"结构 {structure_file.name} 执行异常: {e}")
                    results.append({
                        'structure_name': structure_file.name,
                        'structure_file': str(structure_file),
                        'success': False,
                        'error': str(e)
                    })

        return results

    def _run_single_structure(self, structure_file: Path) -> Dict[str, Any]:
        """
        执行单个结构的Pipeline

        Parameters
        ----------
        structure_file : Path
            结构文件路径

        Returns
        -------
        Dict
            计算结果
        """
        result = {
            'structure_name': structure_file.stem,
            'structure_file': str(structure_file),
            'success': False,
        }

        try:
            # 创建工作目录
            work_dir = self.work_root / structure_file.stem
            if self.pressure_label:
                work_dir = work_dir / self.pressure_label
            work_dir.mkdir(parents=True, exist_ok=True)

            result['work_dir'] = str(work_dir)

            # 创建Pipeline实例
            pipeline = self.pipeline_class(
                structure_file=structure_file,
                work_dir=work_dir,
                **self.pipeline_kwargs
            )

            # 执行Pipeline
            success = pipeline.run()
            result['success'] = success

            # 尝试提取结果
            if success:
                try:
                    from vasp.pipelines.utils import extract_final_energy
                    energy = extract_final_energy(work_dir / "01_relax")
                    if energy:
                        result['energy'] = energy
                except:
                    pass

        except Exception as e:
            logger.error(f"执行Pipeline失败: {e}", exc_info=True)
            result['error'] = str(e)

        return result
