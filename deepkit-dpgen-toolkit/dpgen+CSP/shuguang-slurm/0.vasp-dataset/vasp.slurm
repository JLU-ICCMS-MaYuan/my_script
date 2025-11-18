#!/bin/bash
#SBATCH -J vasp-test    ##作业名
#SBATCH -p tyhcnormal    ##队列
#SBATCH -N 1    ##申请计算节点数
#SBATCH --ntasks-per-node=64  ##每节点进程数

module purge
module load compiler/intel/2017.5.239
module load mpi/intelmpi/2017.4.239

export MKL_DEBUG_CPU_TYPE=5 #加速代码
export MKL_CBWR=AVX2 #使cpu默认支持avx2
export I_MPI_PIN_DOMAIN=numa #内存位置与cpu位置绑定，加速内存读取。对于内存带宽要求高的计算提速明显

for i in 1 2 3 4 5; do
cp INCAR_$i INCAR
srun --mpi=pmi2 /work/home/acvm651ob1/soft/vasp.6.3.2/bin/vasp_std > vasp.log_$i 2>&1 
cp CONTCAR POSCAR
done
