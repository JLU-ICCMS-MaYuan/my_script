"""
文件操作工具

提供统一的文件和目录操作功能。

作者：Claude
创建时间：2025-11-20
"""

import os
import shutil
from pathlib import Path
from typing import List, Optional, Union


def ensure_dir(directory: Union[str, Path]) -> Path:
    """
    确保目录存在，不存在则创建

    Parameters
    ----------
    directory : str or Path
        目录路径

    Returns
    -------
    Path
        目录路径对象
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def copy_file(src: Union[str, Path], dst: Union[str, Path], overwrite: bool = False) -> Path:
    """
    复制文件

    Parameters
    ----------
    src : str or Path
        源文件路径
    dst : str or Path
        目标文件路径
    overwrite : bool
        是否覆盖已存在的文件

    Returns
    -------
    Path
        目标文件路径

    Raises
    ------
    FileExistsError
        如果目标文件已存在且overwrite=False
    """
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        raise FileNotFoundError(f"源文件不存在: {src}")

    if dst.exists() and not overwrite:
        raise FileExistsError(f"目标文件已存在: {dst}")

    # 确保目标目录存在
    dst.parent.mkdir(parents=True, exist_ok=True)

    shutil.copy2(src, dst)
    return dst


def copy_files(src_list: List[Union[str, Path]], dst_dir: Union[str, Path], overwrite: bool = False) -> List[Path]:
    """
    批量复制文件到目标目录

    Parameters
    ----------
    src_list : List[str or Path]
        源文件列表
    dst_dir : str or Path
        目标目录
    overwrite : bool
        是否覆盖已存在的文件

    Returns
    -------
    List[Path]
        目标文件路径列表
    """
    dst_dir = ensure_dir(dst_dir)
    dst_list = []

    for src in src_list:
        src = Path(src)
        dst = dst_dir / src.name
        copy_file(src, dst, overwrite=overwrite)
        dst_list.append(dst)

    return dst_list


def move_file(src: Union[str, Path], dst: Union[str, Path], overwrite: bool = False) -> Path:
    """
    移动文件

    Parameters
    ----------
    src : str or Path
        源文件路径
    dst : str or Path
        目标文件路径
    overwrite : bool
        是否覆盖已存在的文件

    Returns
    -------
    Path
        目标文件路径
    """
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        raise FileNotFoundError(f"源文件不存在: {src}")

    if dst.exists() and not overwrite:
        raise FileExistsError(f"目标文件已存在: {dst}")

    # 确保目标目录存在
    dst.parent.mkdir(parents=True, exist_ok=True)

    shutil.move(str(src), str(dst))
    return dst


def remove_file(file_path: Union[str, Path], ignore_missing: bool = True):
    """
    删除文件

    Parameters
    ----------
    file_path : str or Path
        文件路径
    ignore_missing : bool
        如果文件不存在，是否忽略错误
    """
    file_path = Path(file_path)

    if file_path.exists():
        file_path.unlink()
    elif not ignore_missing:
        raise FileNotFoundError(f"文件不存在: {file_path}")


def remove_dir(dir_path: Union[str, Path], ignore_missing: bool = True):
    """
    删除目录及其内容

    Parameters
    ----------
    dir_path : str or Path
        目录路径
    ignore_missing : bool
        如果目录不存在，是否忽略错误
    """
    dir_path = Path(dir_path)

    if dir_path.exists():
        shutil.rmtree(dir_path)
    elif not ignore_missing:
        raise FileNotFoundError(f"目录不存在: {dir_path}")


def clean_tmp_files(work_dir: Union[str, Path], patterns: Optional[List[str]] = None):
    """
    清理临时文件

    Parameters
    ----------
    work_dir : str or Path
        工作目录
    patterns : List[str], optional
        要清理的文件模式列表（glob格式），默认清理tmp/目录
    """
    work_dir = Path(work_dir)

    if patterns is None:
        # 默认清理tmp目录
        tmp_dir = work_dir / 'tmp'
        if tmp_dir.exists():
            remove_dir(tmp_dir)
    else:
        # 清理指定模式的文件
        for pattern in patterns:
            for file_path in work_dir.glob(pattern):
                if file_path.is_file():
                    remove_file(file_path)
                elif file_path.is_dir():
                    remove_dir(file_path)


def find_files(directory: Union[str, Path], pattern: str, recursive: bool = False) -> List[Path]:
    """
    查找文件

    Parameters
    ----------
    directory : str or Path
        搜索目录
    pattern : str
        文件名模式（glob格式，如 '*.out'）
    recursive : bool
        是否递归搜索子目录

    Returns
    -------
    List[Path]
        找到的文件路径列表
    """
    directory = Path(directory)

    if recursive:
        return list(directory.rglob(pattern))
    else:
        return list(directory.glob(pattern))


def get_latest_file(directory: Union[str, Path], pattern: str) -> Optional[Path]:
    """
    获取最新的文件（按修改时间）

    Parameters
    ----------
    directory : str or Path
        搜索目录
    pattern : str
        文件名模式

    Returns
    -------
    Path or None
        最新文件路径，如果没有找到则返回None
    """
    files = find_files(directory, pattern, recursive=False)

    if not files:
        return None

    return max(files, key=lambda p: p.stat().st_mtime)


def backup_file(file_path: Union[str, Path], backup_suffix: str = '.bak') -> Path:
    """
    备份文件

    Parameters
    ----------
    file_path : str or Path
        要备份的文件
    backup_suffix : str
        备份文件后缀

    Returns
    -------
    Path
        备份文件路径
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"文件不存在: {file_path}")

    backup_path = file_path.with_suffix(file_path.suffix + backup_suffix)

    # 如果备份文件已存在，添加数字后缀
    counter = 1
    while backup_path.exists():
        backup_path = file_path.with_suffix(f"{file_path.suffix}{backup_suffix}.{counter}")
        counter += 1

    shutil.copy2(file_path, backup_path)
    return backup_path
