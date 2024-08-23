#!/bin/sh                           
#------ slurm option --------#
#SBATCH --partition=mpc
#SBATCH --account=hp240139
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=24:00:00

ulimit -s unlimited

#------- Program execution -------#
module load intelmpi/impi_23.2.0
module load intel/23.02.1

conda activate /home/h240012/soft/deepmd-kit

python calypso_run_model_devi.py --all_models ../gen_stru_analy.000/graph.000.pb ../gen_stru_analy.000/graph.001.pb ../gen_stru_analy.000/graph.002.pb ../gen_stru_analy.000/graph.003.pb --type_map Ce Sc H
