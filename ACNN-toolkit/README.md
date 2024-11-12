华为机器检查作业：
```shell
djob
dshowjob
djob -s 504043
djob -T 504043
```

运行任务之前执行：
```shell
conda activate
token_get
```

准备赝势
```shell
python 3-makepot.py
```

准备初始数据集
```shell
conda activate
python /home/share/jildxwlxyljlstdui/home/lijx/playground/torchdemo/interface/calypso2/resource/convert/recycling.py init_data data vasp OUTCAR
```

提交任务命令
```shell
nohup ./scheduling.sh > log 2>&1 &
```

