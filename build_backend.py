import os
import shutil
import sys
import tempfile
from typing import Iterable, Tuple, Union

from setuptools.build_meta import (
    _BuildMetaBackend,
    _file_with_extension,
    no_install_setup_requires,
)


class _PatchedBackend(_BuildMetaBackend):
    """Fallback-friendly backend，避免跨设备重命名导致的安装失败。"""

    def _build_with_temp_dir(
        self,
        setup_command: Iterable[str],
        result_extension: Union[str, Tuple[str, ...]],
        result_directory: str,
        config_settings,
        arbitrary_args: Iterable[str] = (),
    ):
        result_directory = os.path.abspath(result_directory)
        os.makedirs(result_directory, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix=".tmp-", dir=result_directory) as tmp_dist_dir:
            sys.argv = [
                *sys.argv[:1],
                *self._global_args(config_settings),
                *setup_command,
                "--dist-dir",
                tmp_dist_dir,
                *arbitrary_args,
            ]
            with no_install_setup_requires():
                self.run_setup()

            result_basename = _file_with_extension(tmp_dist_dir, result_extension)
            result_path = os.path.join(result_directory, result_basename)
            if os.path.exists(result_path):
                os.remove(result_path)
            # 使用 move 以兼容跨设备/overlay 场景。
            shutil.move(os.path.join(tmp_dist_dir, result_basename), result_path)

        return result_basename


_backend = _PatchedBackend()


def get_requires_for_build_wheel(config_settings=None):
    return _backend.get_requires_for_build_wheel(config_settings)


def get_requires_for_build_editable(config_settings=None):
    return _backend.get_requires_for_build_editable(config_settings)


def prepare_metadata_for_build_wheel(metadata_directory, config_settings=None):
    return _backend.prepare_metadata_for_build_wheel(metadata_directory, config_settings)


def prepare_metadata_for_build_editable(metadata_directory, config_settings=None):
    return _backend.prepare_metadata_for_build_editable(metadata_directory, config_settings)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _backend.build_wheel(wheel_directory, config_settings, metadata_directory)


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    return _backend.build_editable(wheel_directory, config_settings, metadata_directory)
