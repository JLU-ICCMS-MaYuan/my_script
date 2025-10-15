#!/bin/sh                           
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=intel6430
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                         

source /data/home/mayuan/bin/env_gcc-9.2.0

ulimit -s unlimited

cd Base

airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &

airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &

airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &

airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH &
airss.pl -build -max 50000 -seed CeScH 

cd ..