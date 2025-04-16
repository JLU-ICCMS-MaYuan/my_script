#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --no-requeue
#SBATCH --mem=56G
#SBATCH --time=728:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --job-name=train
#SBATCH --output=train-out
source /public/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -v unlimited
export SLURM_EXPORT_ENV=ALL
source /public/home/mayuan/miniconda3/bin/activate /public/home/mayuan/miniconda3/envs/magus

mpirun -n 56 mlp train pot.mtp train.cfg --trained-pot-name=pot.mtp --max-iter=200 --energy-weight=1.0 --force-weight=0.01 --stress-weight=0.001 --scale-by-force=0.0 --weighting=structures --update-mindist --ignore-weights=True
