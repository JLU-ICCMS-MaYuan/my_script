# deepkit+结构预测

##  <span style="color:red">  0. 安装软件，准备目录以及相应的输入文件和训练集

###  <span style="color:yellow"> 0.1 安装deepmd-kit
```shell
https://github.com/deepmodeling/deepmd-kit/releases/tag/v2.2.2
# 离线安装cpu版本
deepmd-kit-2.2.2-cpu-Linux-x86_64.sh
# 下载好之后直接
bash deepmd-kit-2.2.2-cpu-Linux-x86_64.sh
```

###  <span style="color:yellow"> 0.2 安装dpgen
```shell
https://github.com/wangzyphysics/dpgen/tree/devel
# 适用于CALYPSO+DPGEN的版本
# 下载时切记要用git clone 下载，千万不要下载安装包
pip install -e .
```


```shell
  × Getting requirements to build wheel did not run successfully.
  │ exit code: 1
  ╰─> [35 lines of output]
      Traceback (most recent call last):
        File "/work/home/mayuan/miniconda3/envs/deepmd/lib/python3.10/site-packages/pip/_vendor/pyproject_hooks/_in_process/_in_process.py", line 353, in <module>
          main()
        File "/work/home/mayuan/miniconda3/envs/deepmd/lib/python3.10/site-packages/pip/_vendor/pyproject_hooks/_in_process/_in_process.py", line 335, in main
          json_out['return_val'] = hook(**hook_input['kwargs'])
        File "/work/home/mayuan/miniconda3/envs/deepmd/lib/python3.10/site-packages/pip/_vendor/pyproject_hooks/_in_process/_in_process.py", line 118, in get_requires_for_build_wheel
          return hook(config_settings)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/build_meta.py", line 325, in get_requires_for_build_wheel
          return self._get_build_requires(config_settings, requirements=['wheel'])
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/build_meta.py", line 295, in _get_build_requires
          self.run_setup()
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/build_meta.py", line 311, in run_setup
          exec(code, locals())
        File "<string>", line 1, in <module>
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/__init__.py", line 103, in setup
          return distutils.core.setup(**attrs)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/_distutils/core.py", line 147, in setup
          _setup_distribution = dist = klass(attrs)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/dist.py", line 303, in __init__
          _Distribution.__init__(self, dist_attrs)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/_distutils/dist.py", line 283, in __init__
          self.finalize_options()
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools/dist.py", line 654, in finalize_options
          ep(self)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools_scm/_integration/setuptools.py", line 121, in infer_version
          _assign_version(dist, config)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools_scm/_integration/setuptools.py", line 56, in _assign_version
          _version_missing(config)
        File "/tmp/pip-build-env-_j3dooqw/overlay/lib/python3.10/site-packages/setuptools_scm/_get_version_impl.py", line 112, in _version_missing
          raise LookupError(
      LookupError: setuptools-scm was unable to detect version for /work/home/mayuan/mysoftware/dpgen-devel.
      
      Make sure you're either building from a fully intact git repository or PyPI tarballs. Most other sources (such as GitHub's tarballs, a git checkout without the .git folder) don't contain the necessary metadata and will not work.
      
      For example, if you're using pip, instead of https://github.com/user/proj/archive/master.zip use git+https://github.com/user/proj.git#egg=proj
      [end of output]
  
  note: This error originates from a subprocess, and is likely not a problem with pip.
error: subprocess-exited-with-error

```

###  <span style="color:yellow"> deepgen的命令
```shell
nohup dpgen run param.json machine.json 1>log 2>err &
# param.json 参数设置文件
# machine.json 机器配置文件
ps -ef | grep dpgen
# -ef 显示进程具体信息
```


#我修改的部分

```python
/public/home/mayuan/miniconda3/envs/deepmd/lib/python3.10/site-packages/dpdispatcher/machines/pbs.py

#ret, stdin, stdout, stderr = self.context.block_call("qstat -x " + job_id)
ret, stdin, stdout, stderr = self.context.block_call("qstat -l " + job_id)
```