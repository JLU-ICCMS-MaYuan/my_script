(my_scripts) [mayuan@login01 AIRSS]$ cat create.sh 
#!/bin/bash
#SBATCH --job-name=caly_process_unfinished
#SBATCH --partition=intel6430
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1

source /work/home/mayuan/bin/env_gcc.9.5.sh

for i in {1..1000}; do
echo "$i"
./dyn_gcs composition.dat 64
find . -name "*.res" -print0 | xargs -0 -P4 -I {} mv -n {} ../Base
done
