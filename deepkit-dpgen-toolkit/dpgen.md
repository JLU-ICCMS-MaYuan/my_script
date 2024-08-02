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

###  <span style="color:yellow"> 我这里给出一个完美的环境本版
```shell

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