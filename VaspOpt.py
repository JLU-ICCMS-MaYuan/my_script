#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
'''
The function of this scripts is opt many structures at the same time.
The requirement is:
    1. this elements must be made of the same elements, the distinction is only that
       the number of atoms for each element is different.
    2. I have no idea at the moment
The Program Logic:
    1. get all the structure files in the specified directory.
    2. create a new directory with the same name for each structure file
    3. copy the POTCAR of the specified path to every structure directory
    4. Enter each directory to create a Slurm job script
       submit the task of running VASP structure optimization
Use examples:
    首先必须保证目录下的文件是.vasp文件
    VaspOpt.py -d La1B1H10/ -posf -potf -p 200
        -d La1B1H10/          提供结构所在目录，不可以提供所在的文件位置
        -posf                 创建结构优化子目录
        -potf                 是否选择赝势并创建赝势
    
     VaspOpt.py -d all_struct/80.0-90.0/ -inc3 -p 200 -encut 650
        -d La1B1H10/   提供结构所在母目录
        -inc3                创建3步优化的INCARs
        -p 200               设置压强
        -encut 650           设置截断能
    VaspOpt.py -d La1B1H10/ -incf -p 200 -encut 650
        -d La1B1H10/        提供结构所在母目录
        -incf                     创建精细优化的INCAR-fine
        -p 200                    设置压强
        -encut 650                设置截断能
     VaspOpt.py -d all_struct/80.0-90.0/ -rv3
        -d all_struct/80.0-90.0/  提供结构所在母目录
        -rv3                      多步优化
    VaspOpt.py -d La1B1H10/ -rvf
        -d all_struct/80.0-90.0/  提供结构所在母目录
        -rvf                      精细优化
    VaspOpt.py -d 80.0-90.0/ -contf
        -d 80.0-90.0/             提供结构所在母目录
        -contf                    在执行命令的目录下创建一个contcars的目录来存储所有结构优化后的CONTCAR

'''
import os
import re
import shutil
from argparse import ArgumentParser
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

def write_incar1(incar_dirpath, pressure, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_1")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))
        incar.write("PREC     = LOW  \n") 
        incar.write("NCORE    = 4    \n")         
        incar.write("KSPACING = 0.8  \n")            
        incar.write("ISMEAR   = 1    \n")   
        incar.write("SIGMA    = 0.2  \n")   
        incar.write("NELM     = 90   \n")   
        incar.write("NELMIN   = 6    \n")   
        incar.write("EDIFF    = 1e-2 \n")
        incar.write("EDIFFG   = -0.5 \n") 
        incar.write("NSW      = 200  \n")   
        incar.write("IBRION   = 2    \n")   
        incar.write("ISIF     = 2    \n")   
        incar.write("PSTRESS  = {}0  \n".format(str(pressure)))   

def write_incar2(incar_dirpath, pressure, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_2")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))        
        incar.write("PREC     = Normal\n") 
        incar.write("NCORE    = 4    \n")         
        incar.write("KSPACING = 0.5  \n")
        incar.write("ISMEAR   = 1    \n" )   
        incar.write("SIGMA    = 0.2  \n" )   
        incar.write("NELM     = 90   \n" )   
        incar.write("NELMIN   = 6    \n" )   
        incar.write("EDIFF    = 1e-4 \n" )
        incar.write("EDIFFG   = -0.1 \n" )  
        incar.write("NSW      = 200  \n" )   
        incar.write("IBRION   = 2    \n" )   
        incar.write("ISIF     = 4    \n" )   
        incar.write("POTIM    = 0.50 \n" )        
        incar.write("PSTRESS  = {}0  \n".format(str(pressure)) )   

def write_incar3(incar_dirpath, pressure, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_3")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))        
        incar.write("PREC     = A    \n")
        incar.write("NCORE    = 4    \n")         
        incar.write("KSPACING = 0.20 \n")            
        incar.write("ISMEAR   = 1    \n")   
        incar.write("SIGMA    = 0.01 \n")   
        incar.write("NELM     = 90   \n")   
        incar.write("NELMIN   = 6    \n")   
        incar.write("EDIFF    = 1e-8 \n")
        incar.write("EDIFFG   = -0.001\n")
        incar.write("NSW      = 200  \n")   
        incar.write("IBRION   = 2    \n")   
        incar.write("ISIF     = 3    \n")        
        incar.write("PSTRESS  = {}0  \n".format(str(pressure)))    

def fine_incar(incar_dirpath, pressure, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_fine")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))        
        incar.write("PREC     = A    \n")
        incar.write("NCORE    = 4    \n")         
        incar.write("KSPACING = 0.20 \n") 
        incar.write("ISMEAR   = 1    \n")   
        incar.write("SIGMA    = 0.01 \n")   
        incar.write("NELM     = 100  \n")   
        incar.write("NELMIN   = 6    \n")   
        incar.write("EDIFF    = 1e-8 \n")
        incar.write("NSW      = 200  \n")   
        incar.write("IBRION   = 2    \n")   
        incar.write("ISIF     = 3    \n")    
        incar.write("PSTRESS  = {}0  \n".format(str(pressure)))   

def slurm3opt(slurm_dirpath):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurm3opt.sh")
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write("#!/bin/sh                            \n")     
        slurm.write("#SBATCH  --job-name=opt3steps        \n")                         
        slurm.write("#SBATCH  --output=opt3steps.out.%j   \n")                       
        slurm.write("#SBATCH  --error=opt3steps.err.%j    \n")                      
        slurm.write("#SBATCH  --partition=lhy           \n")                   
        slurm.write("#SBATCH  --nodes=1                   \n")             
        slurm.write("#SBATCH  --ntasks=48                 \n")               
        slurm.write("#SBATCH  --ntasks-per-node=48        \n")                        
        slurm.write("#SBATCH  --cpus-per-task=1           \n")                     
        slurm.write("\n\n")
        slurm.write("source /work/env/intel2018           \n")  
        slurm.write("\n\n")
        slurm.write('cp INCAR_1 INCAR                                                 \n')       
        slurm.write('cp POSCAR POSCAR-0                                               \n')         
        slurm.write('mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')                        
        slurm.write('cp -f CONTCAR CONTCAR-1 &&  cp -f CONTCAR POSCAR                 \n')
        slurm.write('echo "opt 1 finished"\n')                                                                          
        slurm.write('cp INCAR_2 INCAR                                                 \n')       
        slurm.write('mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')                        
        slurm.write('cp -f CONTCAR CONTCAR-2 &&  cp -f CONTCAR POSCAR                 \n')                                                                           
        slurm.write('echo "opt 2 finished"\n') 
        slurm.write('cp INCAR_3 INCAR                                                 \n')       
        slurm.write('mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')                        
        slurm.write('cp -f CONTCAR CONTCAR-3 &&  cp -f CONTCAR POSCAR                 \n')                                       
        slurm.write('echo "opt 3 finished"\n')

def slurmFopt(slurm_dirpath):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurmFopt.sh")
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write('#!/bin/sh                                                                \n')     
        slurm.write('#SBATCH  --job-name=opt_fine                                             \n')                         
        slurm.write('#SBATCH  --output=opt_fine.out.%j                                        \n')                       
        slurm.write('#SBATCH  --error=opt_fine.err.%j                                         \n')                      
        slurm.write('#SBATCH  --partition=lhy                                               \n')    # lhy lbt is both ok                
        slurm.write('#SBATCH  --nodes=1                                                       \n')             
        slurm.write('#SBATCH  --ntasks=48                                                     \n')               
        slurm.write('#SBATCH  --ntasks-per-node=48                                            \n')                        
        slurm.write('#SBATCH  --cpus-per-task=1                                               \n')                     
        slurm.write('\n\n                                                                     \n')
        slurm.write('source /work/env/intel2018                                               \n')
        slurm.write('ulimit -s unlimited                                                      \n')
        slurm.write('export I_MPI_ADJUST_REDUCE=3                                             \n')
        slurm.write('export MPIR_CVAR_COLL_ALIAS_CHECK=0                                      \n')
        slurm.write('\n\n                                                                     \n')
        slurm.write('cp INCAR_fine INCAR                                                      \n')                    
        slurm.write('num=0                                                                    \n')                                                                               
        slurm.write('while true;do                                                            \n')              
        slurm.write('        let num+=1                                                       \n')                   
        slurm.write('        echo "run fine vasp opt-$num"                                    \n')                                      
        slurm.write('        mpirun -n 48 /work/software/vasp.6.1.0/vasp_std  > vasp.log 2>&1 \n')                                                                         
        slurm.write('        cp -f CONTCAR CONTCAR-fine &&  cp -f CONTCAR POSCAR              \n')                                                            
        slurm.write("        rows=`sed -n '/F\=/p' OSZICAR | wc -l`                           \n")                                               
        slurm.write('        echo "rows-$rows"                                                \n')                           
        slurm.write('        if [ "$rows" -eq "1" ];then                                      \n')                                    
        slurm.write('                break                                                    \n')                      
        slurm.write('        fi                                                               \n')           
        slurm.write('done                                                                     \n')          
                                 

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--directory-sourcefile",
        default=None,
        type=str,
        required=True,
        dest="directory_sourcefile",
        help="请输入结构所在的目录"
    )
    parser.add_argument(
        "-posf",
        "--create-poscar",
        action="store_true",
        dest="create_poscar",
        default=False,
        help="是否要创建结构优化目录opt_dir, 并且将.vasp结构文件放进opt_dir里面。"
    )
    parser.add_argument(
        "-inc3",
        "--create-incar-for-3-steps",
        action="store_true",
        dest="create_incar_for_3",
        help="是否进行3步结构优化INCAR文件的产生"
    )
    parser.add_argument(
        "-incf",
        "--create-incar-for-fine-step",
        action="store_true",
        dest="create_incar_for_fine",
        help="是否进行精细结构优化INCAR文件的产生"
    )
    parser.add_argument(
        "-p",
        "--pressure",
        default=None,
        type=int,
        dest="pressure",
        help="在哪个压强下进行结构优化"
    )
    parser.add_argument(
        "-encut",
        "--encut",
        default=None,
        type=int,
        dest="encut",
        help="请输入截断能"
    )
    parser.add_argument(
        "-potf",
        "--potcar-flag",
        default=False,
        action="store_true",
        dest="potcar_file_flag",
        help="请输入POTCAR的路径"
    )
    parser.add_argument(
        "-rv3",
        "--run-vasp-structure-relex-3steps",
        default=False,
        action="store_true",
        dest="run_vasp_3steps",
        help="是否要运行多步vasp结构优化"
    )
    parser.add_argument(
        "-rvf",
        "--run-vasp-high-fine",
        default=False,
        action="store_true",
        dest="run_vasp_high_fine",
        help="是否要运行一步vasp精细结构优化"
    )
    parser.add_argument(
        '-contf',
        '--get-contcar-file-flag',
        action='store_true',
        default=False,
        dest='get_contcar_file_flag',
        help='是否将COTCAR'
    )

    args = parser.parse_args()
    directory_sourcefile  = args.directory_sourcefile   
    create_poscar         = args.create_poscar          # 准备POSCAR 建立运行opt的子目录
    create_incar_for_3    = args.create_incar_for_3     # 准备INCAR
    create_incar_for_fine = args.create_incar_for_fine  # 准备INCAR
    pressure              = args.pressure
    encut                 = args.encut
    potcar_file_flag      = args.potcar_file_flag
    run_vasp_3steps       = args.run_vasp_3steps
    run_vasp_high_fine    = args.run_vasp_high_fine
    get_contcar_file_flag = args.get_contcar_file_flag

    # 准备POSCAR 和 POTCAR 建立运行opt的子目录
    if create_poscar and (directory_sourcefile is not None) and (pressure is not None):                                     # 创建结构优化的目录，需要.vasp文件的位置，压强，potcar
        for root, dirs, files in os.walk(directory_sourcefile):

            if len(files)!=1:
                for f in files:
                    if re.search(r"\.vasp$", f): # 找到所有文件中以.vasp结尾的文件名
                        src_filepath = os.path.join(directory_sourcefile, f) # 未被移动到opt_dir中的vasp结构文件路径
                        # dst_dir 这个目录是放置 number-opt 这个最小子目录的上一层母目录
                        # opt_dir 这个目录很重要，运行结构优化的最小子目录
                        dst_dir = os.path.join(directory_sourcefile, f.split(".")[0])
                        opt_dir = os.path.join(dst_dir, str(pressure)+"-opt")
                        pot_lib = os.path.join(dst_dir, "potcar_lib")
                        if not os.path.exists(opt_dir):
                            os.makedirs(opt_dir)
                        if not os.path.exists(pot_lib): # 提前准备好一个赝势库的空目录  
                                os.makedirs(pot_lib)
                        # 在opt_dir目录中准备好POSCAR结构文件
                        poscar_filepath = os.path.join(opt_dir, "POSCAR")      #  被移动到opt_dir中的vasp结构文件经过拷贝后变成POSCAR的文件路径
                        # 创建POSCAR
                        shutil.copy(src_filepath, poscar_filepath)
                        # 并将POSCAR改变为原胞
                        struct  = Structure.from_file(poscar_filepath)
                        spa     = SpacegroupAnalyzer(struct)
                        pstruct = spa.get_primitive_standard_structure()
                        bstruct = spa.get_conventional_standard_structure()
                        Poscar(pstruct).write_file(os.path.join(opt_dir, "PPOSCAR"))
                        Poscar(bstruct).write_file(os.path.join(opt_dir, "BPOSCAR"))
                        shutil.copy(
                            os.path.join(opt_dir, "PPOSCAR"),
                            os.path.join(opt_dir, "POSCAR"),
                        )
                        if potcar_file_flag:
                            # 准备分离的赝势
                            all_element = [ ele.name for ele in pstruct.types_of_species ]
                            POT_dir = os.path.abspath("/work/home/may/POT/vasp_pot1/potpaw_PBE54")
                                                      
                            for element in all_element:
                                print("\n\n ------------element={}-------------".format(element))
                                dirs_list = os.listdir(POT_dir)
                                for d in dirs_list:
                                    if element in d:
                                        print(d)
                                input_valid_flag = False
                                while not input_valid_flag:
                                    potcar_dir = input("请选择你要使用的POTCAR: \n")
                                    if potcar_dir in dirs_list:
                                        src_potcar =os.path.join(POT_dir, potcar_dir, "POTCAR")
                                        dst_potcar = os.path.join(pot_lib, element)
                                        if os.path.exists(src_potcar):
                                            shutil.copy(src_potcar, dst_potcar)
                                            worklog = os.path.join(pot_lib, "pot_selected.log")
                                            with open(worklog, "a") as w:
                                                w.write("{} {}\n".format(potcar_dir, src_potcar))
                                            input_valid_flag = True
                                        else:
                                            input_valid_flag = False
                                    else:
                                        input_valid_flag = False
                            # 合并分立的POTCAR
                            pot_lib = os.path.join(dst_dir, "potcar_lib")
                            pot_name_list = os.listdir(pot_lib)
                            if not os.path.exists(pot_lib):
                                raise FileExistsError("potcar_lib doesn't exist!")
                            for root, dirs, files in os.walk(opt_dir):
                                if "POSCAR" in files:
                                    src_poscar = os.path.join(root, "POSCAR")
                                    elements = open(src_poscar, "r").readlines()[5].split()
                                    print("element {}".format(elements))
                                    src_potcar_path_list = []
                                    for ele in elements:
                                        if ele in pot_name_list:
                                            src_potcar = os.path.join(pot_lib, ele)
                                            if os.path.exists(src_potcar):  
                                                src_potcar_path_list.append(src_potcar)
                                    print("src_potcar_path_list={}".format(src_potcar_path_list))
                                    dst_potcar = os.path.join(opt_dir, "POTCAR")
                                    # 将多个POTCAR写入总的POTCAR中
                                    f = open(dst_potcar, "w")
                                    for potcar_indi in src_potcar_path_list:
                                        for line in open(potcar_indi):
                                            f.writelines(line)
                                    f.close()
                                    os.system("grep TITEL {}".format(dst_potcar))
                                    print("完成了一个POTCAR\n\n\n")
            elif len(files)==1:
                if re.search(r"\.vasp$", files[0]): # 这仅有的一个结构是否以.vasp结尾
                    src_filepath = os.path.join(directory_sourcefile, files[0]) # 未被移动到opt_dir中的vasp结构文件路径
                    # opt_dir 这个目录很重要，运行结构优化的最小子目录
                    opt_dir = os.path.join(directory_sourcefile, str(pressure)+"-opt")
                    if not os.path.exists(opt_dir):
                        os.makedirs(opt_dir)
                    # 在opt_dir目录中准备好.vasp结构文件
                    poscar_filepath = os.path.join(opt_dir, "POSCAR")      #  被移动到opt_dir中的vasp结构文件经过拷贝后变成POSCAR的文件路径
                    # 创建POSCAR
                    shutil.copy(src_filepath, poscar_filepath)
                    # 并将POSCAR改变为原胞
                    struct  = Structure.from_file(poscar_filepath)
                    spa     = SpacegroupAnalyzer(struct)
                    pstruct = spa.get_primitive_standard_structure()
                    bstruct = spa.get_conventional_standard_structure()
                    Poscar(pstruct).write_file(os.path.join(opt_dir, "PPOSCAR"))
                    Poscar(bstruct).write_file(os.path.join(opt_dir, "BPOSCAR"))
                    shutil.copy(
                        os.path.join(opt_dir, "PPOSCAR"),
                        os.path.join(opt_dir, "POSCAR"),
                    )
                    if potcar_file_flag:
                        # 准备分离的赝势
                        all_element = [ ele.name for ele in pstruct.types_of_species ]
                        POT_dir = os.path.abspath("/work/home/may/POT/vasp_pot1/potpaw_PBE54")
                        # 提前准备好一个赝势库的空目录
                        pot_lib = os.path.join(directory_sourcefile, "potcar_lib")
                        if not os.path.exists(pot_lib):
                            os.makedirs(pot_lib)
                        for element in all_element:
                            print("\n\n ------------element={}-------------".format(element))
                            dirs_list = os.listdir(POT_dir)
                            for d in dirs_list:
                                if element in d:
                                    print(d)
                            input_valid_flag = False
                            while not input_valid_flag:
                                potcar_dir = input("请选择你要使用的POTCAR: \n")
                                if potcar_dir in dirs_list:
                                    src_potcar =os.path.join(POT_dir, potcar_dir, "POTCAR")
                                    dst_potcar = os.path.join(pot_lib, element)
                                    if os.path.exists(src_potcar):
                                        shutil.copy(src_potcar, dst_potcar)
                                        worklog = os.path.join(pot_lib, "pot_selected.log")
                                        with open(worklog, "a") as w:
                                            w.write("{} {}\n".format(potcar_dir, src_potcar))
                                        input_valid_flag = True
                                    else:
                                        input_valid_flag = False
                                else:
                                    input_valid_flag = False
                        # 合并分立的POTCAR
                        pot_lib = os.path.join(directory_sourcefile, "potcar_lib")
                        pot_name_list = os.listdir(pot_lib)
                        if not os.path.exists(pot_lib):
                            raise FileExistsError("potcar_lib doesn't exist!")
                        for root, dirs, files in os.walk(opt_dir):
                            if "POSCAR" in files:
                                src_poscar = os.path.join(root, "POSCAR")
                                elements = open(src_poscar, "r").readlines()[5].split()
                                print("element {}".format(elements))
                                src_potcar_path_list = []
                                for ele in elements:
                                    if ele in pot_name_list:
                                        src_potcar = os.path.join(pot_lib, ele)
                                        if os.path.exists(src_potcar):  
                                            src_potcar_path_list.append(src_potcar)
                                print("src_potcar_path_list={}".format(src_potcar_path_list))
                                dst_potcar = os.path.join(opt_dir, "POTCAR")
                                f = open(dst_potcar, "w")
                                for potcar_indi in src_potcar_path_list:
                                    for line in open(potcar_indi):
                                        f.writelines(line)
                                f.close()
                                os.system("grep TITEL {}".format(dst_potcar))
                                print("完成了一个POTCAR")

    #是否进行3步结构优化INCAR文件的产生         
    if create_incar_for_3 and pressure is not None and encut is not None  and directory_sourcefile is not None:
        for root, dirs, files in os.walk(directory_sourcefile):
            if "POSCAR" in files:
                opt_dir = os.path.abspath(root)
                # 创建INCAR
                write_incar1(opt_dir, pressure, encut)
                write_incar2(opt_dir, pressure, encut)
                write_incar3(opt_dir, pressure, encut)

    #是否进行精细结构优化INCAR文件的产生    
    if create_incar_for_fine and pressure is not None and encut is not None and directory_sourcefile is not None:
        for root, dirs, files in os.walk(directory_sourcefile):
            if "POSCAR"  in files: 
                opt_dir = os.path.abspath(root)
                slurmFopt(opt_dir)
                fine_incar(opt_dir, pressure, encut)

    # 是否运行3步结构优化
    if run_vasp_3steps and directory_sourcefile is not None:
        for root, dirs, files in os.walk(directory_sourcefile):
            if "POSCAR"  in files and \
               "POTCAR"  in files and \
               "INCAR_1" in files and \
               "INCAR_2" in files and \
               "INCAR_3" in files :
                opt_dir = os.path.abspath(root)
                # 创建提交作业的脚本 并执行
                slurm3opt(opt_dir)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurm3opt.sh")
                os.chdir(cwd)
    
    # 是否进行1步精细结构优化
    if run_vasp_high_fine and directory_sourcefile is not None:
        for root, dirs, files in os.walk(directory_sourcefile):
            if "POSCAR"     in files and \
               "POTCAR"     in files and \
               "INCAR_fine" in files:
                opt_dir = os.path.abspath(root)
                # 创建提交作业的脚本 并执行
                slurmFopt(opt_dir)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmFopt.sh")
                os.chdir(cwd)
    
    # 是否获得CONTCAR
    if get_contcar_file_flag and directory_sourcefile is not None:
        cwd = os.getcwd(); os.mkdir("contcars"); contcars_dir = os.path.join(cwd, "contcars")
        for root, dirs, files in os.walk(directory_sourcefile):
            if "CONTCAR" in files:
                src_contcar_file = os.path.join(root, "CONTCAR")
                if os.path.exists(src_contcar_file):
                    des_contcar_file = os.path.join(contcars_dir, root.split("/")[-1]+".vasp")
                    shutil.copy(src_contcar_file, des_contcar_file)
                    
                



