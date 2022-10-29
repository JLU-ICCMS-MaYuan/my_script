# qebin_path = "/work/home/may/software/qe-7.1/bin/"
# qe_source_libs = "/work/home/may/POT/qe-pp/all_pbe_UPF_v1.5"
# eliashberg_x_path = "/work/home/may/code/my_script/qe/eliashberg/eliashberg.x"

qebin_path     = "/home/mayuan/mysoftware/q-e-qe-7.1/bin"
qe_source_libs = "/home/mayuan/mysoftware/all_pbe_UPF_v1.5"
eliashberg_x_path = "/work/home/may/code/my_script/qe/eliashberg/eliashberg.x"

# qebin_path = "/public/home/mayuan/software/qe-7.1/bin/"
# qe_source_libs = "/public/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5"
# eliashberg_x_path = "/public/home/mayuan/code/my_script/qe/eliashberg/eliashberg.x"

# qebin_path = "/work/home/mayuan/mysoftware/qe-7.0/bin/"
# intel_compiler = "/work/home/mayuan/intel/oneapi/setvars.sh --force"
# qe_source_libs = "/work/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5"
# eliashberg_x_path = "/work/home/mayuan/code/my_script/qe/eliashberg/eliashberg.x"

bashtitle = '''#!/bin/sh                                                 

source /opt/intel/oneapi/setvars.sh --force      
ulimit -s unlimited            

'''

slurmtitle = '''#!/bin/sh                           
#SBATCH  --job-name=mayuan                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=lhy          
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                          
                          
source /work/home/may/intel/oneapi/setvars.sh --force      
#source /work/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -s unlimited            
                      
'''

pbstitle = '''            
#!/bin/sh                       
#PBS -N    relax                                    
#PBS -q    liuhy         
#PBS -l    nodes=1:ppn=28               
#PBS -j    oe                                      
#PBS -V                                         
                 
source /public/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited        
cd $PBS_O_WORKDIR                  
killall -9 pw.x ph.x

'''