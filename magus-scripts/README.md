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
### 2. 记得修改input.yaml中的pressure，单位GPa
### 3. 记得修改input.yaml中的ppLabel，例如：ppLabel: ['_GW','_sv','']
### 4. 不管是input.yaml的MLCalculator还是MainCalculator, 都要记得设置preProcessing，例如：
```shell
 preProcessing: |            
  source /work/home/may/intel/oneapi/setvars.sh --force
  source activate /work/home/may/miniconda3/envs/magus 
  export I_MPI_ADJUST_REDUCE=3
  export MPIR_CVAR_COLL_ALIAS_CHECK=0
  ulimit -s unlimited
```
### 5. input.yaml中的段落前空格非常重要
### 6. 记得修改input.yaml中的minNAtoms和maxNAtoms，特别是当你做定组分结构预测的时候，一定要保证你的组分的
### 7. 记得激活环境
```shell
source activate /work/home/may/miniconda3/envs/magus
```