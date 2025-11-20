"""
批量绘图和HTML报告生成器

作者：Claude
创建时间：2025-11-19
"""

from pathlib import Path
from typing import List, Dict
import base64
from io import BytesIO
import matplotlib.pyplot as plt


class BatchPlotReporter:
    """
    批量绘图报告生成器

    功能：
    1. 为多个结构批量绘图
    2. 生成HTML汇总报告
    3. 支持图片预览和下载
    """

    def __init__(self, output_dir: Path):
        """
        初始化报告生成器

        Parameters
        ----------
        output_dir : Path
            输出目录
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.plot_data: List[Dict] = []

    def add_plot(self, structure_name: str, plot_type: str, fig: plt.Figure):
        """
        添加一个图表到报告

        Parameters
        ----------
        structure_name : str
            结构名称
        plot_type : str
            图表类型（phonon_band, electron_dos等）
        fig : Figure
            matplotlib图形对象
        """
        # 保存PNG文件
        png_path = self.output_dir / f"{structure_name}_{plot_type}.png"
        fig.savefig(png_path, dpi=300, bbox_inches='tight')

        # 转换为Base64（用于HTML嵌入）
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.read()).decode()

        self.plot_data.append({
            'structure': structure_name,
            'type': plot_type,
            'png_path': png_path,
            'base64': img_base64
        })

        plt.close(fig)

    def generate_html_report(self, output_path: Optional[Path] = None) -> Path:
        """
        生成HTML报告

        Parameters
        ----------
        output_path : Path, optional
            输出HTML文件路径

        Returns
        -------
        Path
            生成的HTML文件路径
        """
        if output_path is None:
            output_path = self.output_dir / 'report.html'

        # 按结构分组
        structures = {}
        for item in self.plot_data:
            struct = item['structure']
            if struct not in structures:
                structures[struct] = []
            structures[struct].append(item)

        # 生成HTML
        html_content = self._generate_html_content(structures)

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"✓ HTML报告已生成: {output_path}")
        return output_path

    def _generate_html_content(self, structures: Dict) -> str:
        """生成HTML内容"""
        html = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>QE计算结果报告</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            text-align: center;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .structure-section {{
            background: white;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .structure-title {{
            font-size: 24px;
            font-weight: bold;
            color: #333;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
        }}
        .plot-item {{
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
            transition: transform 0.2s;
        }}
        .plot-item:hover {{
            transform: scale(1.02);
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        .plot-title {{
            background: #667eea;
            color: white;
            padding: 10px;
            font-weight: bold;
            text-align: center;
        }}
        .plot-image {{
            width: 100%;
            height: auto;
            display: block;
        }}
        .download-btn {{
            display: block;
            text-align: center;
            padding: 8px;
            background: #764ba2;
            color: white;
            text-decoration: none;
            transition: background 0.3s;
        }}
        .download-btn:hover {{
            background: #667eea;
        }}
        .summary {{
            background: #e8f5e9;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 QE计算结果报告</h1>
        <p>批量绘图结果汇总</p>
    </div>

    <div class="summary">
        <h2>📊 统计信息</h2>
        <p>结构数量: {len(structures)}</p>
        <p>图表总数: {len(self.plot_data)}</p>
    </div>
"""

        # 为每个结构生成图表展示
        for struct_name, plots in structures.items():
            html += f"""
    <div class="structure-section">
        <div class="structure-title">📁 {struct_name}</div>
        <div class="plot-grid">
"""
            for plot in plots:
                plot_type_name = self._format_plot_type(plot['type'])
                html += f"""
            <div class="plot-item">
                <div class="plot-title">{plot_type_name}</div>
                <img src="data:image/png;base64,{plot['base64']}" class="plot-image" alt="{plot_type_name}">
                <a href="{plot['png_path'].name}" class="download-btn" download>💾 下载PNG</a>
            </div>
"""
            html += """
        </div>
    </div>
"""

        html += """
</body>
</html>
"""
        return html

    @staticmethod
    def _format_plot_type(plot_type: str) -> str:
        """格式化图表类型名称"""
        type_names = {
            'phonon_band': '声子能带',
            'phonon_dos': '声子态密度',
            'electron_band': '电子能带',
            'electron_dos': '电子态密度',
            'projected_band': '投影能带',
            'projected_dos': '投影态密度',
        }
        return type_names.get(plot_type, plot_type)


# 使用示例
if __name__ == "__main__":
    # 创建报告生成器
    reporter = BatchPlotReporter(output_dir=Path("./reports"))

    # 添加图表（示例）
    # fig1 = plt.figure()
    # ... 绘图代码 ...
    # reporter.add_plot("H3S", "phonon_band", fig1)

    # 生成HTML报告
    # reporter.generate_html_report()
    pass
