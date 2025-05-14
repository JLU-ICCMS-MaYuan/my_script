qebin_path = "/work/home/mayuan/software/qe-7.1/bin"
epwbin_path = "/work/home/mayuan/software/qe-7.1/EPW/bin"
qe_pseudopotential_dir = "/work/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5"
eliashberg_x_path = "/work/home/mayuan/code/my_script/qe/eliashberg/eliashberg.x"

vaspstd_path = "/work/home/mayuan/software/vasp.6.1.0/bin/vasp_std"
vaspgam_path = "/work/home/mayuan/software/vasp.6.1.0/bin/vasp_gam"

potcar_dir = "/work/home/mayuan/POT/vasp_pot1/potpaw_PBE54"


bashtitle = '''#!/bin/sh   
source /opt/intel/oneapi/setvars.sh --force  
ulimit -s unlimited
'''

slurmtitle = '''#!/bin/sh                           
#SBATCH  --job-name=vasp
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
##SBATCH  --partition=intel6240r_384
#SBATCH  --partition=intel6240r_192
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=48                          
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         

##These parameters of slurm is for many nodes parallel.
##You use 4 nodes, occupy 8 cores for each node, the rest cores for restore the data.
##exclusive means only you can use these 4 nodes, nobody can use them with you together.
##SBATCH  --job-name=mayuan
##SBATCH  --output=log.out                       
##SBATCH  --error=log.err                       
##SBATCH  --partition=public
##SBATCH  --nodes=4
##SBATCH  --ntasks=32
##SBATCH  --ntasks-per-node=8
##SBATCH  --cpus-per-task=1                         
##SBATCH  --exclude=node98
##SBATCH  --exclusive

## mpirun -np 32 /work/home/mayuan/software/qe-7.1/EPW/bin/epw.x -npool 32 <epw_elph.in> epw_elph.out

source /work/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
'''

pbstitle = '''#!/bin/sh                   
#PBS -N    mayqe                                    
#PBS -q    liuhy         
#PBS -l    nodes=1:ppn=28               
#PBS -j    oe                                      
#PBS -V  
source /public/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited        
cd $PBS_O_WORKDIR                  
#killall -9 pw.x ph.x
'''

lsftitle = '''#!/bin/bash
#BSUB -n 56
#BSUB -q normal
#BSUB -J myjob
#BSUB -R 'select[ncpus>=56] span[ptile=56]'
#BSUB -o operation.log

##EXAMPLE for many nodes 
##BSUB -n 72
##BSUB -q normal
##BSUB -J myjob
##BSUB -R 'select[ncpus>=56] span[ptile=18]'
##BSUB -o operation.log
## 72/ptile = 72/18 = 4:  4代表你申请的节点数, 18代表你每个节点使用的核数，剩余的核数用来存储


'''

if __name__ == "__main__":
    
    from pathlib import Path

    qebin_path          = Path(qebin_path)
    qe_pseudopotential_dir_path = Path(qe_pseudopotential_dir)
    eliashberg_x_path   = Path(eliashberg_x_path)

    if Path(qebin_path).exists():
        print(f"The path {qebin_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {qebin_path} you set wrong !!!")
    
    if Path(epwbin_path).exists():
        print(f"The path {epwbin_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {epwbin_path} you set wrong !!!")

    if Path(qe_pseudopotential_dir_path).exists():
        print(f"The path {qe_pseudopotential_dir_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {qe_pseudopotential_dir_path} you set wrong !!!")
    if Path(eliashberg_x_path).exists():
        print(f"The path {eliashberg_x_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {eliashberg_x_path} you set wrong !!!")

    vaspstd_path       = Path(vaspstd_path)
    potcar_dir_path    = Path(potcar_dir)
    if Path(vaspstd_path).exists():
        print(f"The path {vaspstd_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {vaspstd_path} you set wrong !!!")
    
    vaspgam_path       = Path(vaspgam_path)
    potcar_dir_path    = Path(potcar_dir)
    if Path(vaspgam_path).exists():
        print(f"The path {vaspgam_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {vaspgam_path} you set wrong !!!")

    if Path(potcar_dir_path).exists():
        print(f"The path {potcar_dir_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {potcar_dir_path} you set wrong !!!")
        
    print("bash      的设置, 请注意:")
    print("    1. 请注意source编译器的设置")
    print(bashtitle, "\n")

    print("slurmtitle 的设置, 请注意:")
    print("    1. 请注意source编译器的设置")
    print("    2. 请注意job-name, partition, ntasks, ntasks-per-node设置")
    print("    3. 请注意是否需要添加以下设置")
    print("        ulimit -s unlimited ")
    print("        export I_MPI_ADJUST_REDUCE=3")
    print("        export MPIR_CVAR_COLL_ALIAS_CHECK=0")
    print(slurmtitle, "\n")

    print("pbstitle 的设置, 请注意:")
    print("    1. 请注意source编译器的设置")
    print("    2. 请注意PBS -N , PBS -q, PBS -l 设置")
    print("    3. 请注意是否需要添加以下设置")
    print("        ulimit -s unlimited")      
    print("        cd $PBS_O_WORKDIR  ")
    print("        export I_MPI_ADJUST_REDUCE=3")
    print("        export MPIR_CVAR_COLL_ALIAS_CHECK=0")
    print(pbstitle, "\n")

    def Write_Tobin(my_scriptrc_ini, qebin_path, vaspbin_path, epwbin_path):
        with open(my_scriptrc_ini, "r") as f:
            content = f.read()
    
        with open(qebin_path, "w") as qe:
            qe.writelines(content)
    
        with open(vaspbin_path, "w") as vasp:
            vasp.write(content)
        
        with open(epwbin_path, "w") as vasp:
            vasp.write(content)
    
    my_scriptrc_ini = Path.home().joinpath(".my_scriptrc.py")
    qebin_path      = Path.home().joinpath("code/my_script/qe/qebin.py")
    vaspbin_path    = Path.home().joinpath("code/my_script/vasp/vaspbin.py")
    epwbin_path     = Path.home().joinpath("code/my_script/epw/epwbin.py")
    
    Write_Tobin(my_scriptrc_ini, qebin_path, vaspbin_path, epwbin_path, )
