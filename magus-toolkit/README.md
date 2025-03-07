## <span style="font-size: 30px; color: lightgreen;"> 网络经验贴
https://www.bilibili.com/read/cv23297817/

我的帖子: https://zhuanlan.zhihu.com/p/685715050


## <span style="font-size: 30px; color: lightgreen;"> 将OUTCAR转化为train.cfg
```shell
mlp convert-cfg OUTCAR train.cfg --input-format=vasp-outcar 
```

## <span style="font-size: 30px; color: lightgreen;"> 检查train.cfg的最小距离
```shell
mlp mindist train.cfg
```

## <span style="font-size: 30px; color: lightgreen;"> 基本流程
### <span style="font-size: 25px; color: red;"> 1. 准备输入文件
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


## <span style="font-size: 30px; color: lightgreen;"> 注意事项
### <span style="font-size: 25px; color: red;"> 1. 记得准备inputFold/MTP/train.cfg

### <span style="font-size: 25px; color: red;"> 2. input.yaml的设置注意事项
<span style="font-size: 18px; color: skyblue;"> **2.1 记得修改input.yaml中的pressure，单位GPa** 

<span style="font-size: 18px; color: skyblue;">  **2.2 记得修改input.yaml中的ppLabel**

例如：ppLabel: ['_GW','_sv','']


<span style="font-size: 18px; color: skyblue;">  **2.3 不管是input.yaml的MLCalculator还是MainCalculator, 都要记得设置preProcessing，**

例如：
```shell
 preProcessing: |            
  source /work/home/may/intel/oneapi/setvars.sh --force
  source activate /work/home/may/miniconda3/envs/magus 
  export I_MPI_ADJUST_REDUCE=3
  export MPIR_CVAR_COLL_ALIAS_CHECK=0
  ulimit -s unlimited
```

<span style="font-size: 18px; color: skyblue;">  **2.4 input.yaml中的段落前空格非常重要**

<span style="font-size: 18px; color: skyblue;">  **2.5 记得修改input.yaml中的minNAtoms和maxNAtoms, 或者min_n_atoms和max_n_atoms**

特别是当你做定组分结构预测的时候，一定要保证你的组分的原子数足够大

<span style="font-size: 18px; color: skyblue;">  **2.6 formula_pool 这个参数非常重要，它控制了产生结构的配比，如果你不主动删除它，它是不会更新的。**

所以它还有一种奇特的用法：你可以手动指定formula_pool中的内容以保证只产生你需要的配比。


<span style="font-size: 18px; color: skyblue;">  **2.7 kill_time的设置

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

<span style="font-size: 18px; color: skyblue;">  **2.8. magus中关于距离的设置**

明确关于原子间距离和体积设置的参数有：
```yaml
d_ratio: 0.7
distance_matrix: [[1.89,1.87,1.35],[1.87,1.85,1.33],[1.35,1.33,0.8]]
radius: [1.37,1.32,0.58]
volume_ratio: 1.1
min_dist: 0.5
```

**摘自magus/utils.py的get_threshold_dict函数和get_distance_dict函数, 其中关于产生结构的代码在magus/generators/random.py**

1. 如果没有设置`radius`，将会根据`radius = [covalent_radii[atomic_numbers[atom]] for atom in symbols]`设置原子半径
2. 如果没有设置`distance_matrix`，将会根据`distance_dict[(si, sj)] = distance_dict[(sj, si)] = d_ratio * (ri + rj)`来设置。其中ri和rj原子半径，并且`threshold_dict`将会按照`d_ratio`设置。
3. 如果设置了`distance_matrix，distance_dict[(si, sj)] = distance_dict[(sj, si)] = distance_matrix[i][j]`，并且`threshold_dict`将会根据`threshold_dict[(si, sj)] = threshold_dict[(sj, si)] = distance_matrix[i][j] / (ri + rj)`设置
4. 关于体积的设置是 volume-ratio of each structure is calculated by cell_volume / SUM(atom_ball_volume). 在定组分中命名为volume_ratio，在变组分中命名为volRatio

<span style="font-size: 18px; color: skyblue;">  **2.9 Dict参数对应的原子间距离**
```shell
grep "Distance Dict"  tem.log
```



### <span style="font-size: 25px; color: red;"> 3. 记得激活环境
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


### <span style="font-size: 25px; color: red;">  5. magus生成结构： 
```
# 读取输入文件生成10个结构
magus generate -i input.yaml -n 10

# 处理结果
magus summary gen.traj

# 保存结构
magus summary gen.traj -s -o poscars

# 指定保存的结构数量
magus summary gen.traj -s -o poscars -n 500
```

### <span style="font-size: 25px; color: red;"> 6. 如何续算magus

一般来说，magus分为`Initialize`, `Generation 1`, `Generation 2`....这么几步。要想续算，必须保证Initialize完成才能续算，不然它永远会从Initialize开始

```shell
============== Initialize ==============
...
============= Generation 1 =============
...
============= Generation 2 =============
...
============= Generation 3 =============
...
```

续算的时候，只需要在原本的命令中加入-r即可。
```shell
nohup magus search -i input.yaml -m -r > tem.log 2>&1 &
echo $! > taskids
```

### <span style="font-size: 25px; color: red;"> 8.magus的卸载
```shell
pip uninstall magus-kit
```

### <span style="font-size: 25px; color: red;"> 9. 在用summary模式时技巧

summary如何执行？`~/magus-master/magus/entrypoints/summary.py`代码位置
```python
class Summary:
    ...
    ...
    def summary(self, filenames, show_number=20, need_sorted=True, sorted_by='Default', reverse=True, save=False, outdir=None):
        filenames = convert_glob(filenames)
        self.prepare_data(filenames)
        show_number = min(len(self.all_frames), show_number)
        self.show_features_table(show_number, reverse, need_sorted, sorted_by)
        if save:
            self.save_atoms(show_number, outdir)
        if self.formula_type == 'var':
            self.plot_phase_diagram()
```


```shell
# 可以手动添加信息，例如额外添加体积信息
magus summary  gen.traj -a volume

# -b 设置边界
# -v 设置是否编组分，只有设置了-v才能画图
# -p 设置被保存结构寻找对称性时的精度
magus summary results/good.traj -sb ehull -a ehull  -b CeH4 CaH2 H2 -v -p 0.00001
```



### <span style="font-size: 25px; color: red;"> 10. 修改了代码中并行计算, 有些机器不支持mpirun命令, 支持srun，直接使用穿行计算 或者srun (calculators/mtp.py的249行)
```python
def calc_grade(self):
    # must have: pot.mtp, train.cfg
    log.info('\tstep 01: calculate grade')
    #exeCmd = f"{self.mtp_runner} -n 1 {self.mtp_exe} calc-grade pot.mtp train.cfg train.cfg "\
    #         "temp.cfg --als-filename=A-state.als"
    exeCmd = f"mlp calc-grade pot.mtp train.cfg train.cfg "\
            "temp.cfg --als-filename=A-state.als"
    exitcode = subprocess.call(exeCmd, shell=True)
    if exitcode != 0:
        raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
```

```python
def relax_with_mtp(self):
    #content = f"{self.mtp_runner} -n {self.num_core} {self.mtp_exe} relax mlip.ini "\
    #          f"--pressure={self.pressure} --cfg-filename=to_relax.cfg "\
    #          f"--force-tolerance={self.force_tolerance} --stress-tolerance={self.stress_tolerance} "\
    #          f"--min-dist={self.min_dist} --log=mtp_relax.log "\
    #          f"--save-relaxed=relaxed.cfg\n"\
    #          f"cat relaxed.cfg?* > relaxed.cfg\n"
    content = f"srun {self.mtp_exe} relax mlip.ini "\
              f"--pressure={self.pressure} --cfg-filename=to_relax.cfg "\
              f"--force-tolerance={self.force_tolerance} --stress-tolerance={self.stress_tolerance} "\
              f"--min-dist={self.min_dist} --log=mtp_relax.log "\
              f"--save-relaxed=relaxed.cfg\n"\
              f"cat relaxed.cfg?* > relaxed.cfg\n"
    if self.mode == 'parallel':
        self.J.sub(content, name='relax', file='relax.sh', out='relax-out', err='relax-err')
        self.J.wait_jobs_done(self.wait_time)
        self.J.clear()
        time.sleep(10)
    elif self.mode == 'serial':
        exitcode = subprocess.call(content, shell=True)
        if exitcode != 0:
            raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
```

```python
def scf_(self, calcPop):
    calc_dir = self.calc_dir
    basedir = '{}/epoch{:02d}'.format(calc_dir, 0)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    shutil.copy("{}/mlip.ini".format(self.input_dir), "{}/pot.mtp".format(basedir))
    shutil.copy("{}/pot.mtp".format(self.ml_dir), "{}/pot.mtp".format(basedir))
    shutil.copy("{}/train.cfg".format(self.ml_dir), "{}/train.cfg".format(basedir))
    dump_cfg(calcPop, "{}/to_scf.cfg".format(basedir), self.symbol_to_type)

    #exeCmd = f"{self.mtp_runner} -n {self.num_core} {self.mtp_exe} calc-efs {basedir}/pot.mtp {basedir}/to_scf.cfg {basedir}/scf_out.cfg"
    exeCmd = f"srun --mpi=pmi2  {self.mtp_exe} calc-efs {basedir}/pot.mtp {basedir}/to_scf.cfg {basedir}/scf_out.cfg"
    #exeCmd = f"mlp calc-efs {0}/pot.mtp {0}/to_scf.cfg {0}/scf_out.cfg".format(basedir)
    exitcode = subprocess.call(exeCmd, shell=True)
    if exitcode != 0:
        raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
    scfpop = load_cfg("{}/scf_out.cfg".format(basedir), self.type_to_symbol)
    for atoms in scfpop:
        enthalpy = (atoms.info['energy'] + self.pressure * atoms.get_volume() * GPa) / len(atoms)
        atoms.info['enthalpy'] = round(enthalpy, 6)
    return scfpop
```

```python
def train(self, epoch=None):
    epoch = epoch or self.n_epoch
    nowpath = os.getcwd()
    os.chdir(self.ml_dir)
    if not self.ignore_weights:
        self.reweighting()
    #content = f"{self.mtp_runner} -n {self.num_core} {self.mtp_exe} train "\
    #          f"pot.mtp train.cfg --trained-pot-name=pot.mtp --max-iter={epoch} "\
    #          f"--energy-weight={self.weights[0]} --force-weight={self.weights[1]} --stress-weight={self.weights[2]} "\
    #          f"--scale-by-force={self.scaled_by_force} "\
    #          f"--weighting=structures "\
    #          f"--update-mindist "\
    #          f"--ignore-weights={self.ignore_weights}\n"
    content = f"srun {self.mtp_exe} train "\
              f"pot.mtp train.cfg --trained-pot-name=pot.mtp --max-iter={epoch} "\
              f"--energy-weight={self.weights[0]} --force-weight={self.weights[1]} --stress-weight={self.weights[2]} "\
              f"--scale-by-force={self.scaled_by_force} "\
              f"--weighting=structures "\
              f"--update-mindist "\
              f"--ignore-weights={self.ignore_weights}\n"

    if self.mode == 'parallel':
        self.J.sub(content, name='train', file='train.sh', out='train-out', err='train-err')
        self.J.wait_jobs_done(self.wait_time)
        self.J.clear()
    elif self.mode == 'serial':
        exitcode = subprocess.call(content, shell=True)
        if exitcode != 0:
            raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
    os.chdir(nowpath)
```

```python
def calc_efs(self, frames):
    if isinstance(frames, Atoms):
        frames = [frames]
    nowpath = os.getcwd()
    os.chdir(self.ml_dir)
    dump_cfg(frames, 'tmp.cfg', self.symbol_to_type)
    #exeCmd = f"{self.mtp_runner} -n {self.num_core} {self.mtp_exe} calc-efs pot.mtp tmp.cfg out.cfg"
    exeCmd = f"srun --mpi=pmi2 {self.mtp_exe} calc-efs pot.mtp tmp.cfg out.cfg"
    exitcode = subprocess.call(exeCmd, shell=True)
    if exitcode != 0:
        raise RuntimeError('MTP calc-efs exited with exit code: {}.'.format(exitcode))
    result = load_cfg('out.cfg', self.type_to_symbol)
    os.remove('tmp.cfg')
    os.remove('out.cfg')
    os.chdir(nowpath)
    return result
```

```python
def get_loss(self, frames):
    nowpath = os.getcwd()
    os.chdir(self.ml_dir)
    dump_cfg(frames, 'tmp.cfg', self.symbol_to_type)
    #exeCmd = f"{self.mtp_runner} -n {self.num_core} {self.mtp_exe} calc-errors pot.mtp tmp.cfg | grep 'Average absolute difference' | awk {{'print $5'}}"
    exeCmd = f"srun {self.mtp_exe} calc-errors pot.mtp tmp.cfg | grep 'Average absolute difference' | awk {{'print $5'}}"
    loss = os.popen(exeCmd).readlines()
    mae_energies, r2_energies = float(loss[1]), 0.
    mae_forces, r2_forces = float(loss[2]), 0.
    mae_stress, r2_stress = float(loss[3]), 0.
    os.remove('tmp.cfg')
    os.chdir(nowpath)
    return mae_energies, r2_energies, mae_forces, r2_forces, mae_stress, r2_stress
```

```python
def select(self, pop):
    nowpath = os.getcwd()
    os.chdir(self.ml_dir)
    dump_cfg(pop, "new.cfg", self.symbol_to_type)
    # content = f"{self.mtp_runner} -n 1 {self.mtp_exe} select-add "\
    #           f"pot.mtp train.cfg new.cfg diff.cfg "\
    #           f"--weighting=structures\n"
    content = f"srun --mpi=pmi2 {self.mtp_exe} select-add "\
              f"pot.mtp train.cfg new.cfg diff.cfg "\
              f"--weighting=structures\n"
    if self.mode == 'parallel':
        self.J.sub(content, name='select', file='select.sh', out='select-out', err='select-err')
        self.J.wait_jobs_done(self.wait_time)
        self.J.clear()
        time.sleep(10)
    elif self.mode == 'serial':
        exitcode = subprocess.call(content, shell=True)
        if exitcode != 0:
            raise RuntimeError('MTP exited with exit code: %d.  ' % exitcode)
    diff_frames = load_cfg("diff.cfg", self.type_to_symbol)
    os.chdir(nowpath)
    if isinstance(pop, Population):
        return pop.__class__(diff_frames)
    return diff_frames
```

### <span style="font-size: 25px; color: red;"> 10. 修改了代码中slurm提交VASP任务后检查任务是否完成的部分

删了一行alldone=False. (parallel/queuemanage.py的274行)

### <span style="font-size: 25px; color: red;"> 11. 种子文件制作

#### 11.1 <span style="font-size: 20px; color: lightyellow;">  magus 是如何读入种子文件的呢? 

在`~/magus-master/magus/search/search.py`文件中
```python
class Magus
    ...
    ...
    def read_seeds(self):
        log.info("Reading Seeds ...")
        seed_frames = read_seeds('{}/POSCARS_{}'.format(self.seed_dir, self.curgen))
        seed_frames.extend(read_seeds('{}/seeds_{}.traj'.format(self.seed_dir, self.curgen)))
        seed_pop = self.Population(seed_frames, 'seed', self.curgen)
        for i in range(len(seed_pop)):
            seed_pop[i].info['gen'] = self.curgen
        return seed_pop
```
#### 11.2 <span style="font-size: 20px; color: lightred;">  magus 什么时候开始读取种子文件，这决定了你怎么使用该程序?
1. 首先你可以在log.txt和tmp.log中grep Reading Seeds, 来判断是否读入了Seeds
2. 目前还不知道什么时候读取种子文件. 目前还不知道在程序运行过程中是否可以读取种子文件

#### 11.3 <span style="font-size: 20px; color: lightred;"> 如何读取种子文件?

1. 创建一个叫做Seeds的目录
2. 然后在其中起名POSCARS_m, 代表在第m代读入种子文件
3. 或者也可以起名seeds_m.traj, 代表在第m代读入种子文件

### <span style="font-size: 25px; color: red;"> 12. pot.mtp 是机器学习的初始势函数

`pot.mtp`这个势函数不能乱选，如果你之前没有训练好的势能，那么就用mlip给的未训练的势函数， 他们存放在这里。三元体系一般选择`20.mtp`
```shell

cd mlip-2-master/untrained_mtps
ls
02.mtp      04.mtp      06.mtp      08.mtp      10.mtp      12.mtp      14.mtp      16.mtp      18.mtp      20.mtp      22.mtp      24.mtp      26.mtp      28.mtp      readme.txt
# 数字代表势函数的level，level越高精度和消耗越高，一般单质给16，化合物给18-20
cp ~/code/mlip-2/untrained_mtps/20.mtp inputFold/MTP/pot.mtp 
# 将某一个拷贝为inputFold/MTP/pot.mtp即可

# 然后依照input.yaml中的mindist，修改inputFold/MTP/pot.mtp中的min_dist和species_count即可，species_count代表元素个数，min_dist是最小的合理距离。
```
min_dist在默认值那里设置了0.5，但是mugs会根据实际体系自动调整. 

MTP中给出的在relax时关于min_dist的含义解释为：`--min-dist=<num>: terminate relaxation if atoms come closer than <num>`，当原子间距离小于0.5Angstrom时停止优化. 你可以通过`mlp help relax`查看。

1. 在`mlFold/MTP/train.sh`中设置了关于自动调整最小距离的参数：
```shell
mpirun -n 48 mlp train pot.mtp train.cfg --trained-pot-name=pot.mtp --max-iter=200 --energy-weight=1.0 --force-weight=0.01 --stress-weight=0.001 --scale-by-force=0.0 --weighting=structures --update-mindist --ignore-weights=True
```
2. 但是比较奇怪的是，在`calcFold/MTP/epoch**/relax.sh`中却没有让其自动调整`min-dist`
```shell
mpirun -n 48 mlp relax mlip.ini --pressure=200 --cfg-filename=to_relax.cfg --force-tolerance=0.001 --stress-tolerance=0.01 --min-dist=0.5 --log=mtp_relax.log --save-relaxed=relaxed.cfg
```

### <span style="font-size: 25px; color: red;"> 13. magus-master/magus/parallel/queuemanage.py中新增加了这样的代码用于slurm系统检测任务运行状态。


```python
wait_command = f"salloc {wait_condition} {self.wait_params} -p {self.queue_name} sleep 10"
```

有些slurm系统需要指定用户名需要加上--account=xxx这个参数
```python
wait_command = f"salloc {wait_condition} {self.wait_params} -p {self.queue_name} --account=hp240139 sleep 10"
```

### <span style="font-size: 25px; color: red;"> 14. 很多时候报错是因为MTP结构优化出现问题。

在做Ce-Sr-H体系的结构预测时，发现它经常连第一代都不能完整跑完。于是我想检查第一代每一次主动学习结构优化Ce-Sr-H的结构优化情况。
```shell
# 这是用来结构话的命令，存储在calcFold/MTP/epoch*/relax.sh中
mpirun -n 48 mlp relax mlip.ini --pressure=200 --cfg-filename=to_relax.cfg --force-tolerance=0.001 --stress-tolerance=0.01 --min-dist=0.5 --log=mtp_relax.log --save-relaxed=relaxed.cfg
# mlip.ini 结构优化设置文件
# --cfg-filename=<str>: Read initial configurations from <str>
# --log=<str>: Write relaxation log to <str>
# --min-dist=<num>: terminate relaxation if atoms come closer than <num>
# --save-relaxed=<str>: Save the relaxed configurations to <str>
```
所以我们只要对比`to_relax.cfg`和`relaxed.cfg`两个文件中结构的数量差距就可以知道`to_relax.cfg`中有多少结构被优化好了。直接进入`calcFold/MTP/`目录后，执行下面的命令即可
```shell
check_epoch_relaxcfg.sh
```
