qebin_path = "/data/home/mayuan/soft/qe-7.0/bin"
qe_source_libs = "/data/home/mayuan/POT/qepp/all_pbe_UPF_v1.5"
eliashberg_x_path = "/data/home/mayuan/code/my_script/qe/eliashberg/eliashberg.x"

vaspbin_path = "/data/home/mayuan/soft/vasp.6.1.0/bin/vasp_std"
potcar_source_libs = "/data/home/mayuan/POT/vasppp/potpaw_PBE54"


bashtitle = '''#!/bin/sh   
source /opt/intel/oneapi/setvars.sh --force  
ulimit -s unlimited
'''

slurmtitle = '''#!/bin/sh                           
#SBATCH  --job-name=mayqe                      
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=intel6430
#SBATCH  --nodes=1                          
#SBATCH  --ntasks=64
#SBATCH  --ntasks-per-node=64
#SBATCH  --cpus-per-task=1                         

source /data/home/mayuan/intel/oneapi/setvars.sh

ulimit -s unlimited
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

if __name__ == "__main__":
    
    from pathlib import Path

    qebin_path          = Path(qebin_path)
    qe_source_libs_path = Path(qe_source_libs)
    eliashberg_x_path   = Path(eliashberg_x_path)

    if Path(qebin_path).exists():
        print(f"The path {qebin_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {qebin_path} you set wrong !!!")
    if Path(qe_source_libs_path).exists():
        print(f"The path {qe_source_libs_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {qe_source_libs_path} you set wrong !!!")
    if Path(eliashberg_x_path).exists():
        print(f"The path {eliashberg_x_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {eliashberg_x_path} you set wrong !!!")

    vaspbin_path       = Path(vaspbin_path)
    potcar_source_libs_path = Path(potcar_source_libs)
    if Path(vaspbin_path).exists():
        print(f"The path {vaspbin_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {vaspbin_path} you set wrong !!!")
    if Path(potcar_source_libs_path).exists():
        print(f"The path {potcar_source_libs_path} you set rightly !!!")
    else:
        raise FileNotFoundError(f"The path {potcar_source_libs_path} you set wrong !!!")
        
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

    def Write_Tobin(qebin_path, vaspbin_path, my_scriptrc_ini):
        with open(my_scriptrc_ini, "r") as f:
            content = f.read()
    
        with open(qebin_path, "w") as qe:
            qe.writelines(content)
    
        with open(vaspbin_path, "w") as vasp:
            vasp.write(content)
    
    my_scriptrc_ini = Path.home().joinpath(".my_scriptrc.py")
    qebin_path      = Path.home().joinpath("code/my_script/qe/qebin.py")
    vaspbin_path    = Path.home().joinpath("code/my_script/vasp/vaspbin.py")
    Write_Tobin(qebin_path, vaspbin_path, my_scriptrc_ini)
