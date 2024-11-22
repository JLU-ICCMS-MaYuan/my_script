# ACNN + lammps 主动学习

##  <span style="font-size: 30px; color: lightgreen;"> 如何在华为机器上使用acnn

```shell
#检查作业
djob
dshowjob

#提交任务
djob -s <任务号>

#杀任务
djob -T <任务号>
```

运行任务之前执行：
```shell
conda activate

# 激活密钥
token_get 
```

##  <span style="font-size: 25px; color: lightgreen;"> 准备ACNN的输入文件

### 输入文件包括:
```shell
resource/
├── calc_energy.sh
├── convert_xsf_vasp.py
├── fp
│   ├── 2-convert_POSCAR_data.py
│   ├── 3-makepot.py
│   ├── INCAR
│   ├── KPOINTS
│   ├── POSCAR
│   ├── POT-Be
│   ├── POTCAR
│   ├── POT-Ce
│   ├── POT-H
│   ├── POT-La
│   ├── POT-Th
│   ├── POT-Y
│   └── sub-fp.sh
├── fun
│   └── fun_lines.sh
├── init_dt
├── inv2
├── iter_scheduling_remote.sh
├── md
│   ├── data.nso
│   ├── in.lmp
│   ├── POSCAR_2x2x2_fromsqs
│   ├── POSCAR_2x2x2_ordered
│   ├── sub-md.sh
│   └── write_lammps.py
├── move_and_rename.sh
├── process_add_folders.sh
├── recycling.py
├── server.sh
├── sub-select.sh
└── sub-train.sh
```

### <span style="font-size: 20px; color: red;"> 1. 准备fp目录下的POTCAR: 千万注意赝势顺序和relabel的POSCAR中元素顺序一致
```shell
# 这个脚本用于根据POSCAR中的元素顺序产生一个POTCAR
python 3-makepot.py
# 你需要在fp目录下生成POTCAR
# 注意要确保POTCAR中元素的顺序与需要重新进行DFT relabel的POSCAR中元素顺序一致
```
**<span style="font-size: 15px; color: yellow;"> 如何确保这一点呢?**

从你的初始数据集中拿出来一个*.xsf文件. 
因为所有需要relabel的结构都需要用convert_xsf_vasp.py将xsf格式转化为POSCAR
所以你只需要拿到一个convert_xsf_vasp.py转化后的POSCAR文件`python convert_xsf_vasp.py <xsf的目录> <存放poscar的目录>`
然后把这个POSCAR和3-makepot.py放到同一个目录下执行`python 3-makepot.py`即可.

或者更简单的, 手动合并POTCAR.

acnn在主动学习的relabel的时候,会拷贝fp目录中的全部内容, 
包括按照你指定顺序合并好的POTCAR(注意千万不要在fp目录中嵌套子目录, 会报错)


### <span style="font-size: 20px; color: red;"> 2. 准备初始数据集

将你的分子动力学 OUTCAR 放到vasp-data目录下, 
如果有多个 OUTCAR, 就做好命名的区分即可, 最好标注清楚压强和超胞大小
例如, OUTCAR_200GPa_2x2x2
然后使用下面的脚本将其转化为xsf文件:
 
```shell
conda activate
python recycling.py init_data data vasp OUTCAR
```

提交任务命令
```shell
nohup ./scheduling.sh > log 2>&1 &
```

