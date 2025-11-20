"""
Pytest配置文件

提供测试fixtures和配置。

作者：Claude
创建时间：2025-11-19
"""

import pytest
from pathlib import Path
import tempfile
import shutil


@pytest.fixture
def test_structure_file():
    """提供测试用的结构文件路径"""
    return Path("/home/mayuan/code/my_script/test/H3S/std_H6S2_Im-3m_229_.vasp")


@pytest.fixture
def temp_work_dir():
    """提供临时工作目录"""
    temp_dir = Path(tempfile.mkdtemp(prefix="qe_test_"))
    yield temp_dir
    # 测试后清理
    if temp_dir.exists():
        shutil.rmtree(temp_dir)


@pytest.fixture
def sample_config():
    """提供示例配置"""
    return {
        'ecutwfc': 80,
        'kpoints': '6 6 6',
        'conv_thr': '1.0d-8'
    }
