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


##  <span style="color:yellow"> deepgen+calypso的使用注意事项

###  <span style="color:green"> 关于record.dpgen的解释
```shell
# 第0代
00 检查准备好的输入文件 
01 train
02 analysis
03 motivation and exploration

# 第1代
10 检查准备好的输入文件 
11 train
12 analysis
13 motivation and exploration

...
...

```

如果你要从02开始续算，那么就把02后面的03 10 11 12 13 ... 都删掉了. 这表示重新做motivation and exploration

###  <span style="color:green"> 注意切换dpgen的分支

下载好dpgen的代码后，用这个切换分支
```shell
git checkout origin/devel
```

###  <span style="color:green"> 注意 param.json 一对互斥的参数

```json
# 参数1
"calypso_input_path":"./calypso_input",        
"model_devi_max_iter": 15,
```

```json
# 参数2
"sys_configs": []
```

这两个参数都是可以用来控制dpgen采样的。参数1优先被读取，如果参数1没有设置，就读取参数2. 

至于参数1和参数2有什么不同，参见：https://github.com/wangzyphysics/dpgen/blob/devel/examples/run/dp-calypso-vasp/param.json

###  <span style="color:green"> 注意 calypso_run_model_devi.py 的运行. 

`iter.00000/01.model_devi/model_devi_results/`里面有个脚本`calypso_run_model_devi.py` 只能在主节点上运行, 具体是 在`iter.000000/01.model_devi/model_devi_results/` 目录下运行下面的命令:
```shell
python calypso_run_model_devi.py --all_models ../gen_stru_analy.000/graph.000.pb ../gen_stru_analy.000/graph.001.pb ../gen_stru_analy.000/graph.002.pb ../gen_stru_analy.000/graph.003.pb --type_map Ce Sc H
```

在运行过程中会产生task.000.000这个文件. 

在运行完毕后会在iter.000000/01.model_devi/model_devi_results/中产生一个Model_Devi.out的文件. 

并且在运行完毕后dpgen会自动在iter.000000/01.model_devi/record.calypso中的文件末尾会添加数字4

```shell
1 0
2 
3
4
4
```

然而在riken的主节点机器上运行这个脚本会导致一堆core-python文件的爆出. 所以最好还是把它放到远程机器上运行. 在执行上述`python calypso_run_model_devi.py ....`命令前, 先在iter.000000/01.model_devi/record.calypso中添加数字4, 然后手动执行命令, 然后执行` nohup dpgen run param.json machine.json 1>log 2>err &
`