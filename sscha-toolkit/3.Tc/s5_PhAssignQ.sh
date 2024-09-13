#!/bin/sh                           
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
killall -9 vasp_std; killall -9 pw.x; killall -9 ph.x
echo "run scf.fit"                                                
srun --mpi=pmi2 /work/home/acvm651ob1/soft/qe7_install/bin/pw.x -npool 4 <scffit.in> scffit.out
echo "run scf"                                                    
srun --mpi=pmi2 /work/home/acvm651ob1/soft/qe7_install/bin/pw.x -npool 4 <scf.in> scf.out
echo "run split_ph"                                                    
srun --mpi=pmi2 /work/home/acvm651ob1/soft/qe7_install/bin/ph.x -npool 4 <split_ph.in> split_ph.out