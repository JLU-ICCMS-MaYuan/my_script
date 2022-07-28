#!/public/home/mayuan/miniconda3/envs/cage/bin/python3
'''
Use example:
    VaspPhono.py -d phonopy/ -pposf -scell 2 2 2 
        -d phonopy/                     计算声子的目录路径
        -pposf                          找到原胞
        -scell 2 2 2                    扩包 2 2 2   
    VaspPhono.py -d scripts_tests/phonopy/ -inc-DFPT -encut 400 -k 50 50 50
        -d scripts_tests/phonopy/       计算声子的目录路径
        -inc-DFPT                       创建DFPT计算的INCAR
        -encut 400                      
    VaspPhono.py -d phonopy/ -DFPT
        -d phonopy/                     计算声子的目录路径
        -DFPT                           提交DFPT的计算任务

    VaspPhono.py -d scripts_tests/phonopy/ -inc-Disp -encut 400
        -d scripts_tests/phonopy/       计算声子的目录路径
        -inc-Disp                       创建INCAR
        -encut 400                      设置INCAR截断能

    VaspPhono.py -d phonopy-Disp/ -Disp
        -d phonopy-Disp/                计算声子的目录路径
         -Disp                          提交有限位移的计算任务

    VaspPhono.py -d 600-DFPT/ -pposf -scell 4 4 4 -encut 650 -k 30 30 30 -inc-DFPT
        -d 600-DFPT/ 
        -pposf -scell 4 4 4 
        -encut 650 
        -k 30 30 30 
        -inc-DFPT

    VaspPhono.py -d 500-DFPT/ -DFPT-prog
    
    VaspPhono.py -d 500-Disp/ -pposf -scell 5 5 5 -encut 650 -k 40 40 40 -inc-Disp -Disp
    VaspPhono.py -d 500-Disp/ -Disp-prog
'''

import os
import re
import shutil
import seekpath
from argparse import ArgumentParser
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
# from pymatgen.io.phonopy import get_Displaced_structures

def DFPT_incar(incar_dirpath, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_DFPT")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0      \n")   
        incar.write("ICHARG   = 2      \n")   
        incar.write("ENCUT    = {}     \n".format(str(encut)))        
        incar.write("PREC     = A      \n")
        incar.write("#NCORE    = 4     \n")         
        incar.write("ISMEAR   = 1      \n")   
        incar.write("SIGMA    = 0.05   \n")   
        incar.write("NELM     = 100    \n")   
        incar.write("NELMIN   = 6      \n")   
        incar.write("EDIFF    = 1e-8   \n")
        incar.write("IBRION   = 8      \n")   
        incar.write("IALGO    = 38     \n")
        incar.write("ADDGRID  = .TRUE. \n")
        incar.write("POTIM    = 0.01   \n" ) 
        incar.write("LWAVE    = .FALSE.\n" )  
        incar.write("LCHARG   = .FALSE.\n" ) 
        incar.write("LREAL    = .FALSE.\n" ) 

def Displacemnet_incar(incar_dirpath, encut):
    incar_filepath = os.path.join(incar_dirpath, "INCAR_Displacement")
    with open(incar_filepath, "w") as incar:
        incar.write("ISTART   = 0    \n")   
        incar.write("ICHARG   = 2    \n")   
        incar.write("ENCUT    = {}   \n".format(str(encut)))        
        incar.write("PREC     = A    \n")
        incar.write("#NCORE    = 4   \n")         
        incar.write("ISMEAR   = 0    \n")   
        incar.write("SIGMA    = 0.01 \n")   
        incar.write("NELM     = 100  \n")   
        incar.write("EDIFF    = 1e-8 \n")
        incar.write("IBRION   = -1   \n")   
        incar.write("IALGO    = 38   \n")

        incar.write("LREAL    =.FALSE.\n")
        incar.write("LWAVE    =.FALSE.\n")
        incar.write("LCHARG   =.FALSE.\n")
        incar.write("ADDGRID  = .TRUE.\n")


def slurmDFPT(slurm_dirpath):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurmDFPT.sh")
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write("#!/bin/sh                           \n")     
        slurm.write("#SBATCH  --job-name=ls              \n")                         
        slurm.write("#SBATCH  --output=DFPT.out.%j       \n")                       
        slurm.write("#SBATCH  --error=DFPT.err.%j        \n")                      
        slurm.write("#SBATCH  --partition=lhy          \n")    # lhy lbt is both ok                
        slurm.write("#SBATCH  --nodes=1                  \n")             
        slurm.write("#SBATCH  --ntasks=48                \n")               
        slurm.write("#SBATCH  --ntasks-per-node=48       \n")                        
        slurm.write("#SBATCH  --cpus-per-task=1          \n")                     
        slurm.write("                                  \n\n")
        slurm.write("source /work/env/intel2018          \n")  
        slurm.write("ulimit -s unlimited                 \n")
        slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
        slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
        slurm.write("                                  \n\n")
        slurm.write('echo "run fine DFPT"                                              \n')
        slurm.write('cp -f INCAR_DFPT INCAR                                            \n')
        slurm.write('cp -f SPOSCAR POSCAR                                              \n')
        slurm.write('mpirun -np 48 /work/software/vasp.6.1.0/vasp_std > vasp.log 2>&1  \n')                        

def slurmDisplacement(slurm_dirpath, number):
    slurm_script_filepath = os.path.join(slurm_dirpath, "slurmDisp.sh")
    print(slurm_script_filepath)
    with open(slurm_script_filepath, "w") as slurm:
        slurm.write("#!/bin/sh                           \n")     
        slurm.write("#SBATCH  --job-name=Disp            \n")                         
        slurm.write("#SBATCH  --output=Disp.out.%j       \n")                       
        slurm.write("#SBATCH  --error=Disp.err.%j        \n")                      
        slurm.write("#SBATCH  --partition=lhy          \n")    # lhy lbt is both ok                
        slurm.write("#SBATCH  --nodes=1                  \n")             
        slurm.write("#SBATCH  --ntasks=48                \n")               
        slurm.write("#SBATCH  --ntasks-per-node=48       \n")                        
        slurm.write("#SBATCH  --cpus-per-task=1          \n")                     
        slurm.write("                                  \n\n")
        slurm.write("source /work/env/intel2018          \n")
        slurm.write("ulimit -s unlimited                 \n")
        slurm.write("export I_MPI_ADJUST_REDUCE=3        \n")
        slurm.write("export MPIR_CVAR_COLL_ALIAS_CHECK=0 \n")
        slurm.write("                                  \n\n")
        slurm.write('echo "run Displacement-{}"          \n'.format(number))
        slurm.write('mpirun -np 48 /work/software/vasp.6.1.0/vasp_std > vasp.log 2>&1  \n')                        

def write_DFPT_band_conf(band_conf_dirpath, species, dim,  path_lable_list, path_coords):
    
    from itertools import chain 
    
    band_conf_filepath = os.path.join(band_conf_dirpath, "band.conf")
    with open(band_conf_filepath, "w") as f:
        f.write("FORCE_CONSTANTS=READ    \n")
        f.write("ATOM_NAME={}            \n".format(' '.join(species)))
        f.write("DIM={}                  \n".format(' '.join(dim)))
        f.write("NPOINTS=101             \n")
        f.write("EIGENVECTORS=.TRUE.     \n")
        f.write("BAND_LABELS={}          \n".format(' '.join(path_lable_list)))
        path_coords = list(chain.from_iterable(path_coords)); path_coords=list(map(str, path_coords))
        f.write("BAND={}                 \n".format(' '.join(path_coords)))
                                                             
def write_Disp_band_conf(band_conf_dirpath, species, dim,  path_lable_list, path_coords):
    
    from itertools import chain 
    
    band_conf_filepath = os.path.join(band_conf_dirpath, "band.conf")
    with open(band_conf_filepath, "w") as f:
        f.write("ATOM_NAME={}            \n".format(' '.join(species)))
        f.write("DIM={}                  \n".format(' '.join(dim)))
        f.write("NPOINTS=101             \n")
        f.write("EIGENVECTORS=.TRUE.     \n")
        f.write("BAND_LABELS={}          \n".format(' '.join(path_lable_list)))
        path_coords = list(chain.from_iterable(path_coords)); path_coords=list(map(str, path_coords))
        f.write("BAND={}                 \n".format(' '.join(path_coords)))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        '-d',
        '--vasp-phono-input-directory',
        dest="vasp_phono_input_directory",
        default=None,
        action='store',
        type=str
    )
    parser.add_argument(
        '-pposf',
        '--create-privimite-poscar',
        dest="create_privimite_poscar",
        default=False,
        action='store_true',
        help="是否创建原胞, 并覆盖之前的POSCAR"
    )
    parser.add_argument(
        "-inc-DFPT",
        "--create-incar-for-phono-DFPT",
        action="store_true",
        dest="create_incar_for_phono_DFPT",
        help="是否创建DFPT声子计算INCAR"
    )
    parser.add_argument(
        "-inc-Disp",
        "--create-incar-for-phono-Displacement",
        action="store_true",
        dest="create_incar_for_phono_Displacement",
        help="是否创建有限位移声子计算INCAR"
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
        '-k',
        '--create-KPOINTS',
        dest="create_KPOINTS",
        type=float,
        nargs="+",
        action="store",
        help="产生KPOINTS, 保证每个方向k点数量在2~3个"
    )
    parser.add_argument(
        "-scell",
        "--create-supercell-and-Displaced_structures",
        action='store',
        dest="create_supercell_and_Displaced_structures",
        help="为声子计算生成一组对称不相等的位移结构",
        type=int,
        nargs="+"
    )
    parser.add_argument(
        "-DFPT",
        "--run-vasp-DFPT-phonopy",
        action='store_true',
        dest="run_vasp_DFPT_phonopy",
        default=False,
        help="进行DFPT声子计算",
    )
    parser.add_argument(
        "-Disp",
        "--run-vasp-Displacement-phonopy",
        action='store_true',
        dest="run_vasp_Displacement_phonopy",
        default=False,
        help="进行Displacement声子计算",
    )
    parser.add_argument(
        '-DFPT-prog',
        '--DFPT-progressing-data',
        action='store_true',
        default=False,
        dest='progressing_data_DFPT',
        help="是否处理DFPT计算的结果"
    )
    parser.add_argument(
        '-Disp-prog',
        '--Disp-progressing-data',
        action='store_true',
        default=False,
        dest='progressing_data_Disp',
        help="是否处理Displacement计算的结果"
    )
    args = parser.parse_args()

    vasp_phono_input_directory                = args.vasp_phono_input_directory
    create_privimite_poscar                   = args.create_privimite_poscar
    create_incar_for_phono_DFPT               = args.create_incar_for_phono_DFPT
    create_incar_for_phono_Displacement       = args.create_incar_for_phono_Displacement
    encut                                     = args.encut
    create_KPOINTS                            = args.create_KPOINTS
    create_supercell_and_Displaced_structures = args.create_supercell_and_Displaced_structures
    run_vasp_DFPT_phonopy                     = args.run_vasp_DFPT_phonopy
    run_vasp_Displacement_phonopy             = args.run_vasp_Displacement_phonopy
    progressing_data_DFPT                     = args.progressing_data_DFPT
    progressing_data_Disp                     = args.progressing_data_Disp

    # 准备原胞 PPOSCAR 
    if create_privimite_poscar and (vasp_phono_input_directory is not None):
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files and "POTCAR" in files:
                poscar_file = os.path.join(root, "POSCAR")
                struct  = Structure.from_file(poscar_file)
                pstruct = SpacegroupAnalyzer(struct).get_primitive_standard_structure()
                Poscar(pstruct).write_file(poscar_file)
                
    # 扩包，超胞 SPOSCAR 为声子计算生成一组对称不相等的位移结构
    if create_supercell_and_Displaced_structures is not None:
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files:
                # PLAN A
                cwd = os.getcwd()
                os.chdir(root)
                os.system('phonopy -d --dim="{} {} {}"'.format(
                    create_supercell_and_Displaced_structures[0],
                    create_supercell_and_Displaced_structures[1],
                    create_supercell_and_Displaced_structures[2],
                ))
                os.chdir(cwd)
                with open(os.path.join(root, "work.log"), "w") as f:
                    f.write("{} {} {}".format(
                    create_supercell_and_Displaced_structures[0],
                    create_supercell_and_Displaced_structures[1],
                    create_supercell_and_Displaced_structures[2],
                    ))
                poscar_file = os.path.join(root, "POSCAR")
                poscar_init = os.path.join(root, "POSCAR-init")
                shutil.copy(poscar_file, poscar_init)
                # PLAN B 失败的原因在于 : 不知道该怎么获得phonopy_Disp.yaml文件
                sposcar_file = os.path.join(root, "SPOSCAR")
                sstruct      = Structure.from_file(sposcar_file)
                print("supercell lattice\n", sstruct.lattice)
                # supercell_struct_list = get_Displaced_structures(
                #     struct,
                #     atom_Disp=0.01,
                #     supercell_matrix=create_supercell_and_Displaced_structures
                # )
                # sposcar = os.path.join(root, "SPOSCAR")
                # Poscar(supercell_struct_list[0]).write_file(sposcar)
                # for i in range(1, len(supercell_struct_list)):
                #     number = str(i).rjust(3,'0')  # 方法会返回一个原字符串右对齐,并使用空格填充至长度 width 的新字符串。如果指定的长度小于字符串的长度则返回原字符串
                #     poscar_Disp = os.path.join(root, "POSCAR-"+number)
                #     Poscar(supercell_struct_list[i]).write_file(poscar_Disp)

    # 是否创建DFPT声子计算INCAR
    if create_incar_for_phono_DFPT and (encut is not None) and (vasp_phono_input_directory is not None):
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files and "POTCAR" in files:
                DFPT_incar(root, encut)
            else:
                print("缺少文件")              
    # 是否创建有限位移声子计算INCAR
    if create_incar_for_phono_Displacement and (encut is not None) and (vasp_phono_input_directory is not None):
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files and "POTCAR" in files:
                Displacemnet_incar(root, encut)
            else:
                print("缺少文件")    

    # 是否创建KPOINTS
    if create_KPOINTS is not None:
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files and "POTCAR" in files:
                
                from pymatgen.io.vasp import Kpoints
                sstruct = Structure.from_file(os.path.join(root, "SPOSCAR"))
                kpoints = Kpoints.automatic_density_by_lengths(
                    sstruct, 
                    length_densities=create_KPOINTS,
                    force_gamma=True)
                kpoints.write_file(os.path.join(root, "KPOINTS"))
                print(kpoints)
            else:
                print("缺少文件")              
    
    # 是否进行DFPT声子计算
    if run_vasp_DFPT_phonopy and (vasp_phono_input_directory is not None): 
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if "POSCAR" in files and \
              "SPOSCAR" in files and \
           "INCAR_DFPT" in files and \
              "KPOINTS" in files and \
               "POTCAR" in files :
                # 创建提交作业的脚本 并执行
                shutil.copy(
                    os.path.join(root, "SPOSCAR"),
                    os.path.join(root, "POSCAR")
                )
                slurmDFPT(root)
                cwd = os.getcwd()
                os.chdir(root)
                os.system("sbatch slurmDFPT.sh")
                os.chdir(cwd)  
            else:
                print("缺少文件")              
                
    # 进行Displacement声子计算
    if run_vasp_Displacement_phonopy and (vasp_phono_input_directory is not None): 
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if  "POSCAR"  in files and \
                "SPOSCAR" in files and \
                "INCAR_Displacement" in files and \
                "KPOINTS" in files and \
                "POTCAR" in files :
                # 创建提交作业的脚本 并执行
                patter = re.compile(r"POSCAR\-[0-9]{3}")
                poscar_number_list = [patter.match(x).group() for x in files if patter.match(x)]
                for poscar_number in poscar_number_list:
                        dst_number_dir = os.path.join(root, "disp-" + poscar_number.split("-")[-1])
                        if not os.path.exists(dst_number_dir):
                            os.makedirs(dst_number_dir)
                        src_poscar = os.path.join(root, poscar_number)       ; dst_poscar = os.path.join(dst_number_dir, "POSCAR");   shutil.copy(src_poscar, dst_poscar)
                        src_potcar = os.path.join(root, "POTCAR")            ; dst_potcar = os.path.join(dst_number_dir, "POTCAR");   shutil.copy(src_potcar, dst_potcar)
                        src_incar  = os.path.join(root, "INCAR_Displacement"); dst_incar  = os.path.join(dst_number_dir, "INCAR" );   shutil.copy(src_incar, dst_incar )
                        src_kpoints= os.path.join(root, "KPOINTS")           ; dst_kpoints= os.path.join(dst_number_dir,"KPOINTS");   shutil.copy(src_kpoints, dst_kpoints)
                        slurmDisplacement(dst_number_dir, poscar_number.split("-")[-1])
                        
                        cwd = os.getcwd()
                        os.chdir(dst_number_dir)
                        os.system("sbatch slurmDisp.sh")
                        os.chdir(cwd)
            else:
                print("缺少文件")              
    
    # 处理DFPT数据
    if progressing_data_DFPT and (vasp_phono_input_directory is not None): 
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if  "POSCAR-init"        in files and \
                "phonopy_disp.yaml"  in files and \
                "INCAR_DFPT"         in files and \
                "work.log"           in files and \
                "vasprun.xml"        in files:
                
                cwd = os.getcwd()
                os.chdir(root)
                os.system("phonopy --fc vasprun.xml")
                os.chdir(cwd)

                poscar_file = os.path.join(root, "POSCAR-init")
                struct_pymatgen = Structure.from_file(poscar_file)
                lattice_type    = SpacegroupAnalyzer(struct_pymatgen).get_lattice_type()
                crystal_system  = SpacegroupAnalyzer(struct_pymatgen).get_crystal_system()
                composition     = struct_pymatgen.composition.get_el_amt_dict()
                species         = composition.keys()
                cell            = struct_pymatgen.lattice.matrix
                coords          = struct_pymatgen.frac_coords
                numbers         = struct_pymatgen.atomic_numbers
                dim             = open(os.path.join(root, "work.log"), "r").read()
 
                # 方案一 使用seekpath寻找高对称点
                # structure=(cell, coords, numbers)
                # output=seekpath.get_path(structure,  symprec=1e-04, angle_tolerance=1e-2)
                # print(output.get("path")); input()
                # high_coords = output.get("point_coords")
                # print(high_coords); input()

                # 方案二 使用ase寻找高对称点
                from ase.io import read
                struct_ase      = read(poscar_file)
                lat             = struct_ase.cell.get_bravais_lattice()
                path_lable      = lat.special_path
                path_lable_list = [[ p for p in pp] for pp in path_lable.split(",") ];  special_points = lat.get_special_points()
                path_coords     = [list(special_points[k]) for k in path_lable_list[0]]
                # print(path_lable_list[0], "\n\n", path_coords)
                write_DFPT_band_conf(root, species, dim, path_lable_list[0], path_coords)

                cwd = os.getcwd()
                os.chdir(root)
                os.system("phonopy -p -s band.conf -c POSCAR-init")
                os.system("phonopy-bandplot  --gnuplot> band.dat")
                os.chdir(cwd)
            else:
                print("缺少文件")              

    # 处理Displacement数据
    if progressing_data_Disp and (vasp_phono_input_directory is not None): 
        for root, dirs, files in os.walk(vasp_phono_input_directory):
            if  "POSCAR-init"        in files and \
                "phonopy_disp.yaml"  in files and \
                "INCAR_Displacement" in files and \
                "work.log"           in files and \
                dirs                 is not None  :
                m_num = str(len(dirs)).rjust(3,'0')
                cwd = os.getcwd()
                os.chdir(root)
                os.system("phonopy -f disp-{001..%s}/vasprun.xml" %(m_num))
                os.chdir(cwd)

                poscar_file = os.path.join(root, "POSCAR-init")
                struct_pymatgen = Structure.from_file(poscar_file)
                lattice_type    = SpacegroupAnalyzer(struct_pymatgen).get_lattice_type()
                crystal_system  = SpacegroupAnalyzer(struct_pymatgen).get_crystal_system()
                composition     = struct_pymatgen.composition.get_el_amt_dict()
                species         = composition.keys()
                cell            = struct_pymatgen.lattice.matrix
                coords          = struct_pymatgen.frac_coords
                numbers         = struct_pymatgen.atomic_numbers
                dim             = open(os.path.join(root, "work.log"), "r").read()
 
                # 方案一 使用seekpath寻找高对称点
                # structure=(cell, coords, numbers)
                # output=seekpath.get_path(structure,  symprec=1e-04, angle_tolerance=1e-2)
                # print(output.get("path")); input()
                # high_coords = output.get("point_coords")
                # print(high_coords); input()

                # 方案二 使用ase寻找高对称点
                from ase.io import read
                struct_ase      = read(poscar_file)
                lat             = struct_ase.cell.get_bravais_lattice()
                path_lable      = lat.special_path
                path_lable_list = [[ p for p in pp] for pp in path_lable.split(",") ];  special_points = lat.get_special_points()
                path_coords     = [list(special_points[k]) for k in path_lable_list[0]]
                # print(path_lable_list[0], "\n\n", path_coords)
                write_Disp_band_conf(root, species, dim, path_lable_list[0], path_coords)

                cwd = os.getcwd()
                os.chdir(root)
                os.system("phonopy -p -s band.conf -c POSCAR-init")
                os.system("phonopy-bandplot  --gnuplot> band.dat")
                os.chdir(cwd)
            else:
                print("缺少文件")              
