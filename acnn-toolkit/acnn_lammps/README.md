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

##  <span style="font-size: 25px; color: lightgreen;"> 提交任务
```shell
nohup ./scheduling.sh > log 2>&1 &
```

##  <span style="font-size: 25px; color: lightgreen;"> 测试模型精度

### <span style="font-size: 20px; color: red;"> 1. 进入某一代主动学习
```shell
$ ls
00  01  02  03  04  05  06  07  08 ...  24  25   log  order  resource  scheduling.sh  taskids

$ cd 25

$ ls
2-train-log   Ap_inv_1.bin   Ap_inv_57.bin  in.acnn                    model-restart  select_label  sub-eval.sh      train_dt
3-md          Ap_inv_39.bin  Ap_inv_58.bin  inv2                       plot.png       slct          sub-select.sh
4-select-log  Ap_inv_4.bin   Ap_inv_90.bin  eval-model-10000-train_dt  labeling       pre-select    sub-train.sh
```
### <span style="font-size: 20px; color: red;"> 2. 准备机器学习输入文件
```shell
cp inv2 in.acnn
```

### <span style="font-size: 20px; color: red;"> 3. 提交任务 

将`emodel model-restart/model-10000  train_dt> ss 2>&1`命令写入你的提交任务的脚本`sub-eval.sh`中, 例如：
```shell
#!/bin/sh
#===========================================================
# Configure DSUB resources
#===========================================================
#DSUB --job_type cosched
#DSUB -n MMM
#DSUB -A root.jildxwlxyljlstdui
#DSUB -q root.default
#DSUB -R cpu=16
#DSUB -N 1
##DSUB -pn cccs-share-agent-[117,126,129]
#DSUB -oo donau.out.%J
#DSUB -eo donau.err.%J
#===========================================================
#/home/share/jildxwlxyljlstdui/home/lijx/playground/torchdemo-dev/build/acnn -train ./inv2 > train-log 2>&1
export PATH=/home/share/jildxwlxyljlstdui/home/lijx/software/torchdemo/build2:$PATH

module purge
export OMP_NUM_THREADS=16

module use /home/HPCBase/workspace/public/software/modules
module load compilers/gcc/gcc9.3.0
module load compilers/bisheng/bisheng2.1.0
module load mpi/hmpi/1.1.1/hmpi1.1.1-bisheng2.1.0

ldd $(which acnn)


emodel model-restart/model-10000  train_dt> ss 2>&1
```

然后直接用命令，例如slurm系统中，提交即可。 **注意你可能需要修改`model-restart/model-10000`的模型的名称**
```shell
sbatch sub-eval.sh
```

说明：`emodel`是一个脚本，其内容为如下，里面执行了一个叫`evalan.py`的python脚本, 它可以画图，然后给出数据精度的分布。
```shell
#!/bin/bash
set -e

if [ $# -eq 2 ];then
    MODEL=$1
    DT=$2
else
    echo "usage: emodel <MODEL PATH> <DT PATH>"
fi



sed "s|evalmodelpath *= *.*|evalmodelpath   =   $MODEL|" in.acnn > in.acnn.tmp
sed -i "s|evaldatapath *= *.*|evaldatapath    =   $DT|" in.acnn.tmp


BMODEL=$(basename $MODEL)
BDT=$(basename $DT)

EVAL_FILE="eval-$BMODEL-$BDT"
acnn -eval in.acnn.tmp > $EVAL_FILE

python $(which evalan.py) $EVAL_FILE | sort -gk 4

rm in.acnn.tmp
```
如果最终没有画出图，可以手动执行绘图：`python $(which evalan.py) $EVAL_FILE | sort -gk 4`


### <span style="font-size: 20px; color: red;"> 4. 分析数据集，从中筛选出坏结构

在resource目录下执行下面的命令， 它会将原子收力大于指定值的结构文件从`prefix.xsf`改名为`prefix`
```shell
for i in *_add init_dt; do cd $i; check_dt 50; cd ..; done
```

有时候这种方法还是没有筛选出全部的坏结构。
我们需要进一步分析`emodel model-restart/model-10000  train_dt> ss 2>&1`执行后在ss文件中给出的信息。

比如：`ss`文件中会对所有进行`eval`的数据集按照力排序。我们需要找到其中受力特别大的结构并删除。
以此类推，我们也可以对能量排序，将能量偏差特别大的结构也删除。
