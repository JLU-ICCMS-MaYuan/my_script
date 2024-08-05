## 网络经验贴
1. https://www.bilibili.com/read/cv23297817/
2. 


## 将OUTCAR转化为train.cfg
```shell
mlp convert-cfg OUTCAR train.cfg --input-format=vasp-outcar 
```

## 检查train.cfg的最小距离
```shell
mlp mindist train.cfg
```

## 基本流程
### 1. 准备输入文件
```shell
inputFold
├── MTP
│   ├── mlip.ini
│   └── train.cfg
├── Vasp1
│   └── INCAR
├── Vasp2
│   └── INCAR
└── Vasp3
    └── INCAR

input.yaml
slurmjob
```

```shell
#!/bin/sh                           
#SBATCH  --job-name=magus
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=lhy          
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         
 
source /work/home/may/intel/oneapi/setvars.sh --force
source activate /work/home/may/miniconda3/envs/magus 

magus search -m
```


## 注意事项
### 1. 记得准备inputFold/MTP/train.cfg
### 2. input.yaml的设置注意事项
<span style="font-size: 15px; color: skyblue;"> **2.1 记得修改input.yaml中的pressure，单位GPa** 

<span style="font-size: 15px; color: skyblue;">  **2.2 记得修改input.yaml中的ppLabel**

例如：ppLabel: ['_GW','_sv','']


<span style="font-size: 15px; color: skyblue;">  **2.3 不管是input.yaml的MLCalculator还是MainCalculator, 都要记得设置preProcessing，**

例如：
```shell
 preProcessing: |            
  source /work/home/may/intel/oneapi/setvars.sh --force
  source activate /work/home/may/miniconda3/envs/magus 
  export I_MPI_ADJUST_REDUCE=3
  export MPIR_CVAR_COLL_ALIAS_CHECK=0
  ulimit -s unlimited
```

<span style="font-size: 15px; color: skyblue;">  **2.4 input.yaml中的段落前空格非常重要**

<span style="font-size: 15px; color: skyblue;">  **2.5 记得修改input.yaml中的minNAtoms和maxNAtoms, 或者min_n_atoms和max_n_atoms**

特别是当你做定组分结构预测的时候，一定要保证你的组分的原子数足够大

<span style="font-size: 15px; color: skyblue;">  **2.6 formula_pool 这个参数非常重要，它控制了产生结构的配比，如果你不主动删除它，它是不会更新的。**

所以它还有一种奇特的用法：你可以手动指定formula_pool中的内容以保证只产生你需要的配比。


<span style="font-size: 15px; color: skyblue;">  kill_time的设置

kill_time: 86400 设置好了之后就可以控制slurm，pbs系统提交的作业可以运行的最长时间。

它的代码是在 parallel/queuemanage.py 中
```python
hours = self.kill_time // 3600
minites = (self.kill_time % 3600) // 60
seconds = int(self.kill_time % 60)
...
f'#SBATCH --time={hours}:{minites}:{seconds}\n'
...
```


### 3. 记得激活环境
```shell
source activate /work/home/may/miniconda3/envs/magus
```
这是因为:如果你的magus是在主节点上运行, 那么有一部分MPT的程序要在主节点上运行. 需要激活环境. 具体来说是这一段代码要在主节点运行:
```python
# magus-master/magus/calculators/mtp.py
246     def calc_grade(self):
247         # must have: pot.mtp, train.cfg
248         log.info('\tstep 01: calculate grade')
249         exeCmd = f"{self.mtp_runner} -n 1 {self.mtp_exe} calc-grade pot.mtp train.cfg train.cfg "\
250                  "temp.cfg --als-filename=A-state.als"
251         #exeCmd = f"mlp calc-grade pot.mtp train.cfg train.cfg "\
252         #         "temp.cfg --als-filename=A-state.als"
253         exitcode = subprocess.call(exeCmd, shell=True)
254         if exitcode != 0:
255             raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
# 翻译过来就是运行:
mpirun -n 1 mlp calc-grade pot.mtp train.cfg train.cfg temp.cfg --als-filename=A-state.als

```

### 4. 记得检查输出的文件中的Distance Dict参数对应的原子间距离
```shell
grep "Distance Dict"  tem.log
```
### 5. volume_ratio在magus中的含义是：In our program, volume-ratio of each structure is calculated by cell_volume / SUM(atom_ball_volume). 在定组分中命名为volume_ratio，在变组分中命名为volRatio

### 6. magus生成结构： 
```
# 读取输入文件生成10个结构
magus generate -i input.yaml -n 10

# 处理结果
magus summary gen.traj

# 保存结构
magus summary gen.traj -s -o poscars
```
### 7.magus的卸载
```shell
pip uninstall magus-kit
```

### 8. 在用summary模式时，可以手动添加信息，例如额外添加体积信息：
```shell
magus summary  gen.traj -a volume
```

### 9. 修改了代码中计算grade的部分，不使用mpirun计算，直接使用穿行计算  (calculators/mtp.py的249行)

### 10. 修改了代码中slurm提交VASP任务后检查任务是否完成的部分，删了一行alldone=False. (parallel/queuemanage.py的274行)

### 11. 种子文件制作：创建一个叫做Seeds的文件，然后在其中起名POSCARS_m, 代表在第m代读入种子文件。

### 12. pot.mtp 是机器学习的初始势函数，这个势函数不能乱选，如果你之前没有训练好的势能，那么就用mlip给的未训练的势函数， 他们存放在这里：
```shell

cd mlip-2-master/untrained_mtps
ls
02.mtp      04.mtp      06.mtp      08.mtp      10.mtp      12.mtp      14.mtp      16.mtp      18.mtp      20.mtp      22.mtp      24.mtp      26.mtp      28.mtp      readme.txt
# 数字代表势函数的level，level越高精度和消耗越高，一般单质给16，化合物给18-20
cp ~/code/mlip-2/untrained_mtps/20.mtp inputFold/MTP/pot.mtp 
# 将某一个拷贝为inputFold/MTP/pot.mtp即可

# 然后依照input.yaml中的mindist，修改inputFold/MTP/pot.mtp中的min_dist和species_count即可，species_count代表元素个数，min_dist是最小的合理距离。
```

### 13. magus-master/magus/parallel/queuemanage.py中新增加了这样的代码用于slurm系统检测任务运行状态。


```python
wait_command = f"salloc {wait_condition} {self.wait_params} -p {self.queue_name} sleep 10"
```

有些slurm系统需要指定用户名需要加上--account=xxx这个参数
```python
wait_command = f"salloc {wait_condition} {self.wait_params} -p {self.queue_name} --account=hp240139 sleep 10"
```
