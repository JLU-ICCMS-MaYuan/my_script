"""
Lambda计算文件备份工具

用于备份不同μ*值的lambda计算输出文件。

重构自: mytoolkit/qe_toolkit/backup_lambda.py

作者：Claude
创建时间：2025-11-20
"""

import shutil
from pathlib import Path
from typing import List, Optional

from rich.console import Console


class LambdaBackup:
    """Lambda计算文件备份器"""

    # 需要备份的文件列表
    FILES_TO_BACKUP = [
        'lambda.in',
        'lambda.out',
        'alpha2F.dat',
        'INPUT',
        'ALPHA2F.OUT',
        'ELIASHBERG.OUT',
        'ELIASHBERG_IA.OUT',
        'ELIASHBERG_GAP_T.OUT',
    ]

    def __init__(self, work_dir: Path = None):
        """
        初始化备份器

        Parameters
        ----------
        work_dir : Path, optional
            工作目录，默认当前目录
        """
        self.work_dir = Path(work_dir) if work_dir else Path.cwd()
        self.console = Console()

    def backup_for_mu(self, mu: float) -> List[str]:
        """
        为特定μ*值备份文件

        Parameters
        ----------
        mu : float
            屏蔽常数μ*值

        Returns
        -------
        List[str]
            成功备份的文件列表
        """
        mu_str = f"{mu:.2f}"
        backed_up = []

        self.console.print(f"\n[bold]备份μ* = {mu_str}的lambda计算文件:[/]")

        for filename in self.FILES_TO_BACKUP:
            src_file = self.work_dir / filename
            dst_file = self.work_dir / f"{mu_str}-{filename}"

            if src_file.exists():
                try:
                    shutil.copy2(src_file, dst_file)
                    backed_up.append(filename)
                    self.console.print(f"  ✓ {filename} -> {dst_file.name}", style="green")
                except Exception as e:
                    self.console.print(f"  ✗ {filename} 备份失败: {e}", style="red")
            else:
                self.console.print(f"  - {filename} 不存在", style="dim")

        return backed_up

    def backup_multiple_mu(self, mu_values: List[float]):
        """
        为多个μ*值备份文件

        Parameters
        ----------
        mu_values : List[float]
            μ*值列表
        """
        self.console.print(f"\n[bold cyan]开始备份{len(mu_values)}个μ*值的文件...[/]")

        all_backed_up = {}

        for mu in mu_values:
            backed_up = self.backup_for_mu(mu)
            all_backed_up[mu] = backed_up

        # 打印摘要
        self.console.print("\n[bold]备份摘要:[/]")
        for mu, files in all_backed_up.items():
            self.console.print(f"  μ* = {mu:.2f}: {len(files)}/{len(self.FILES_TO_BACKUP)} 个文件")

    def list_backups(self) -> dict:
        """
        列出所有备份文件

        Returns
        -------
        dict
            {μ*值: [备份文件列表]}
        """
        backups = {}

        for file_path in self.work_dir.iterdir():
            if not file_path.is_file():
                continue

            name = file_path.name

            # 检查是否是备份文件（格式：0.10-lambda.out）
            parts = name.split('-', 1)
            if len(parts) == 2:
                try:
                    mu = float(parts[0])
                    original_name = parts[1]

                    if original_name in self.FILES_TO_BACKUP:
                        if mu not in backups:
                            backups[mu] = []
                        backups[mu].append(original_name)
                except ValueError:
                    continue

        return backups

    def print_backups(self):
        """打印所有备份文件"""
        backups = self.list_backups()

        if not backups:
            self.console.print("\n[yellow]未找到备份文件[/]")
            return

        self.console.print("\n[bold]已有的备份:[/]")

        for mu in sorted(backups.keys()):
            files = backups[mu]
            self.console.print(f"\n[cyan]μ* = {mu:.2f}[/] ({len(files)} 个文件)")

            for filename in sorted(files):
                self.console.print(f"  • {filename}", style="dim")


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(description="备份Lambda计算文件")
    parser.add_argument("mu", type=float, nargs='*',
                        help="μ*值（可指定多个）")
    parser.add_argument("-l", "--list", action='store_true',
                        help="列出已有备份")
    parser.add_argument("-d", "--directory", type=str, default=None,
                        help="工作目录（默认：当前目录）")

    args = parser.parse_args()

    # 创建备份器
    backup = LambdaBackup(work_dir=args.directory)

    # 列出备份
    if args.list:
        backup.print_backups()
        return

    # 执行备份
    if not args.mu:
        backup.console.print("[red]错误: 请指定至少一个μ*值[/]")
        backup.console.print("用法: python lambda_backup.py 0.10 0.13")
        return

    backup.backup_multiple_mu(args.mu)


if __name__ == "__main__":
    main()
