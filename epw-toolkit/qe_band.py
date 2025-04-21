import argparse
import os
import shutil
import sys
import subprocess
import re
import math
import numpy as np
import xml.etree.ElementTree as ET  
from qe_opt import extract_qe_input
from qe_phonon import convert_vasp_kpath_to_qe, qe2vasp


def generate_band_input(data, total_electrons, num_kpoints):
    print("Generate band.in ...")
    qe_input_content=[]
    for block_name in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
        qe_input_content.append(f"&{block_name}")
        for key, value in data[block_name].items():
            if key == 'calculation':
                qe_input_content.append(f"  {key} = 'bands'")
            else:
                qe_input_content.append(f"  {key} = {value}")
        if block_name == "SYSTEM":
            qe_input_content.append(f"  nbnd = {math.ceil(total_electrons/2*1.5) + 10}")
        qe_input_content.append("/")

    qe_input_content.append("ATOMIC_SPECIES")
    for species in data['ATOMIC_SPECIES']:
        element = species['element']
        mass = species['mass']
        pseudopotential = species['pseudopotential']
        qe_input_content.append(f"  {element} {mass} {pseudopotential}")

    cell_parameters = data['CELL_PARAMETERS']
    qe_input_content.append(f"CELL_PARAMETERS {cell_parameters['units']}")
    for vector in cell_parameters['vectors']:
        qe_input_content.append("  "+ " ".join(vector))

    atomic_positions = data['ATOMIC_POSITIONS']
    qe_input_content.append(f"ATOMIC_POSITIONS {atomic_positions['units']}")
    for position in atomic_positions['positions']:
        element = position['element']
        coordinates = " ".join(position['coordinates'])
        qe_input_content.append(f"  {element} {coordinates}")
    
    qe_input_text = "\n".join(qe_input_content)
    with open("band.in",'w') as file:
        file.write(qe_input_text)
        file.write("\n")
        kpath_content = convert_vasp_kpath_to_qe('KPATH.in',num_kpoints)
        
        file.write(f"K_POINTS crystal_b\n")
        file.write("\n".join(kpath_content))
    
    print("Done.")


def extract_electrons_from_pseudopotential(pseudopotential_file):
    z_valence_pattern = re.compile(r'z_valence="\s*([\d.]+)"')

    with open(pseudopotential_file, 'r') as file:
        for line in file:
            match = z_valence_pattern.search(line)
            if match:
                return float(match.group(1))
    raise ValueError(f"z_valence not found in pseudopotential file : {pseudopotential_file}")


def calculate_total_electrons(data, pseudopotentials_dir):
    total_electrons = 0
    for species in data['ATOMIC_SPECIES']:
        element = species['element']
        pseudopotential_file = os.path.join(pseudopotentials_dir, species['pseudopotential'])
        electrons_per_atom = extract_electrons_from_pseudopotential(pseudopotential_file)
        atom_count = sum(1 for position in data['ATOMIC_POSITIONS']['positions'] if position['element'] == element)
        total_electrons += electrons_per_atom * atom_count
    return total_electrons


def generate_bands_input(prefix):
    print("Generate bands.in file ...")
    with open('bands.in','w') as file:
        file.write(f"&BANDS\n")
        file.write(f"  prefix = '{prefix}'\n")
        file.write(f"  lsym = .false.\n")
        file.write(f"/\n")
    print("Done.")


def generate_projwfc_input(prefix):
    print("Generate projwfc.fat.in file ...")
    with open('projwfc.fat.in','w') as file:
        file.write(f"&projwfc\n")
        file.write(f"  outdir = './' \n")
        file.write(f"  prefix = '{prefix}'\n")
        file.write(f"  lsym = .false.\n")
        file.write(f"  filproj = 'fatband'\n")
        file.write(f"/\n")
    print("Done.")


def generate_dos_input(prefix):
    print("Generate dos.in file ...")
    with open('dos.in','w') as file:
        file.write(f"&dos\n")
        file.write(f"  prefix = '{prefix}'\n")
        file.write(f"  outdir = './'\n")
        file.write(f"  bz_sum = 'tetrahedra_opt'\n")
        file.write(f"/\n")
    print("Done.")


def plot_dos(prefix):
     # 初始化列表来存储能量值和DOS值  
    energies = []  
    dos_values = []  
    efermi = None
  
    import matplotlib.pyplot as plt

    # 读取文件并提取能量和DOS数据  
    with open(f'{prefix}.dos', 'r') as file:  
        for line in file:  
            if line.strip() and not line.startswith('#'):  
                parts = line.split()  
                energy = float(parts[0]) - efermi  # 校准能量值  
                dos = float(parts[1])  
                energies.append(energy)  
                dos_values.append(dos)  
            elif line.startswith('#'):
                parts = line.split()
                efermi = float(parts[8])
            else:
                pass

    # 绘制DOS图  
    plt.plot(energies, dos_values)  
    plt.xlabel('Energy - EFermi (eV)')  
    plt.ylabel('DOS')  
    plt.title('Density of States')  
    
    plt.ylim([-0.01,5])
    
    plt.axvline(x=0, color='red', linestyle='--')  # 在EFermi处画一条垂直线  
    plt.savefig('qe_dos.png',dpi=300)
    plt.show()  




def extract_label(file_name,num_kpoints):
    label = []
    nkpoints = []
    index = []
    with open(file_name, 'r') as vasp_file:
        lines = vasp_file.readlines()
    
    for line in lines:
        if line.strip() and "Line-Mode" not in line and "Reciprocal" not in line and "K-Path" not in line:
            parts = line.split()
            if len(parts) == 4:
                index.append(parts[3])
    
    count = 1
    for i in range(1,len(index)):
        if index[i] != index[i-1]:
            count += 1

    for i in range(0, len(index), 2):
        if i < len(index) -2:
            start_k = index[i]
            end_k = index[i+1]
            next_k = index[i+2]
            if start_k =='GAMMA':
                label.append(r"$\Gamma$")
                nkpoints.append(num_kpoints)
            else:
                label.append(start_k)
                nkpoints.append(num_kpoints)
            if next_k == end_k:
                pass
            else:
                if end_k =='GAMMA':
                    label.append(r"$\Gamma$")
                    nkpoints.append(1)
                else:
                    label.append(end_k)
                    nkpoints.append(1)
        else:
            start_k = index[i]
            end_k = index[i+1]
            if start_k == 'GAMMA':
                label.append(r"$\Gamma$")
                nkpoints.append(num_kpoints)
            else:
                label.append(start_k)
                nkpoints.append(num_kpoints)
            if end_k == 'GAMMA':
                label.append(r"$\Gamma$")
                nkpoints.append(1)
            else:
                label.append(end_k)
                nkpoints.append(1)
    
    return label,nkpoints


#extract fermi energy from {prefix}.save file
def extract_fermi(prefix):
    hartree2ev=27.211407953
    file_name = os.path.join(f"{prefix}.save","data-file-schema.xml")
    tree = ET.parse(file_name)
    root = tree.getroot()  
  
    # extract fermi energy if exists
    try:
        fermi_energy = float(root.find('.//fermi_energy').text)*hartree2ev
    except:
        fermi_energy = None
    
    # extract eigenvalues and occupations  
    eigenvalues_list = root.findall('.//eigenvalues')  
    occupations_list = root.findall('.//occupations[@size]')  
  
    eigen_list = []
    occup_list = []
    for eigenvalues, occupations in zip(eigenvalues_list, occupations_list):  
        eigen_list.append(eigenvalues.text.strip().split())  
        occup_list.append(occupations.text.strip().split())

    # judge the system is metal or semiconductor
    band_index = None
    for i in range(len(occup_list[0])):
        if float(occup_list[0][i]) > 0.99 and float(occup_list[0][i+1]) < 0.01:
            band_index = i
            break
    
    if not band_index:
        #this is a metal
        if not fermi_energy:
            #no use smearing method
            print("The system is metal. Don't use fix or tetrahedron occupation.")
            sys.exit()
    
    #check whether exist band gap
    hasgap = True
    for occup in occup_list:
        if float(occup[i]) < 0.99 or float(occup[i+1]) > 0.01:
            hasgap = False
    
    if hasgap:
        #the system is semiconductor, find the valence band bottom
        band_energy = 9999.99
        for eigenvalue in eigen_list:
            band_energy = min(band_energy,float(eigenvalue[band_index+1]))

        fermi_energy = band_energy * hartree2ev
        system = 'semiconductor'
    
    else:
        #the system is metal, check whether find fermi energy
        system = 'metal'
        if not fermi_energy:
            #no use smearing method
            print("The system is metal. Don't use fix or tetrahedron occupation.")
            sys.exit()
    
    return system, fermi_energy
    

#extract band energy from band.out file
def extract_band_energy(file_name):
    hartree2ev=27.211407953
    file_name = os.path.join(f"{prefix}.save","data-file-schema.xml")
    tree = ET.parse(file_name)
    root = tree.getroot()  
  
    # extract eigenvalues
    eigenvalues_list = root.findall('.//eigenvalues')  
  
    eigen_list = []
    for eigenvalues in eigenvalues_list:  
        eigen_list.append(eigenvalues.text.strip().split())  
    
    for i in range(len(eigen_list)):
        for j in range(len(eigen_list[0])):
            eigen_list[i][j] = float(eigen_list[i][j]) * hartree2ev

    return eigen_list 


#extract projection bands
def extract_orbital(file_name):
    atom_index = []
    with open(file_name,'r') as file:
        ischar = False
        next(file)
        for line in file:
            if line.strip().split()[1].isalpha():
                ischar = True
                atom_index.append(line.strip().split()[1])
            elif not any(char.isalpha() for char in line) and ischar:
                break
        
        for line in file:
            if len(line.strip().split()) == 3:
                parts = line.strip().split()
                break
        
        num_proj = int(parts[0])
        num_kpoints = int(parts[1])
        num_bands = int(parts[2])
        next(file)

        projwfc = np.zeros((num_proj,num_kpoints,num_bands))
        index = None
        symbol = None
        proj = None
        index_next = None
        symbol_next = None
        proj_next = None
        index_list = []
        symbol_list = []
        proj_list = []
        for i in range(num_proj):
            line = file.readline()
            parts = line.strip().split()
            if not index and not symbol and not proj:
                index = parts[1]
                symbol = parts[2]
                proj = parts[3]
            else:
                index_next = parts[1]
                symbol_next = parts[2]
                proj_next = parts[3]
            if not index_next and not symbol_next and not proj_next:
                initial = 0
                index_list.append(index)
                symbol_list.append(symbol)
                proj_list.append(proj)
                for j in range(num_kpoints):
                    for k in range(num_bands):
                        projwfc[initial][j][k] = float(file.readline().strip().split()[2])
            elif index==index_next and symbol == symbol_next and proj == proj_next:
                for j in range(num_kpoints):
                    for k in range(num_bands):
                        projwfc[initial][j][k] += float(file.readline().strip().split()[2])
            else:
                initial += 1
                index = index_next
                symbol = symbol_next
                proj = proj_next
                index_list.append(index)
                symbol_list.append(symbol)
                proj_list.append(proj)
                for j in range(num_kpoints):
                    for k in range(num_bands):
                        projwfc[initial][j][k] += float(file.readline().strip().split()[2])
        
    return index_list, symbol_list, proj_list, projwfc


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description='Generate qe band input file.')  
    parser.add_argument('--num_kpoints', type=int, default = 50, help='number of kpoints between high symmetry path.')
    parser.add_argument('--pw_cmd', type=str, default='mpirun -np 20 pw.x -nk 4',help='parallel run command for pw.x run.')

    args = parser.parse_args()  

    # if os.path.exists("qe_opt"):
    #     if os.path.exists("qe_opt/scf.in"):
    #         pass
    #     else:
    #         print("Error! Can not find file qe_opt/scf.in. Please check whether qe_opt.py run successfully.")
    #         sys.exit()
    # else:
    #     print("Error! Can not find qe_opt/ directory. Please check whether qe_opt.py run successfully.")
    #     sys.exit()
    
    # if os.path.exists("band"):
    #     print("band/ directory has exists. Don't do anything.")
    #     sys.exit()

    # print("Create directory band/.")
    # os.mkdir("band")

    #find all upf file 
    # upf_files = [f for f in os.listdir("./") if f.endswith('.upf')]
    # for file_name in upf_files:
    #     shutil.copy(file_name,'./band')
    # shutil.copy("qe_opt/scf.in","./band")

    os.chdir("band")
    #extract qe data
    shutil.copy("../phonon/scf.in","./")   
    qe_input=extract_qe_input("scf.in")

    #generate POSCAR which will used in VASPKIT to generate high symmetry path
    qe2vasp(qe_input)

    #generate high symmetry path
    vaspkit_commands = '3\n303\n'
    subprocess.run(['vaspkit'], input=vaspkit_commands, text=True)

    prefix = qe_input['CONTROL']['prefix'].strip("'")
    atomic_positions = []
    for position in qe_input['ATOMIC_POSITIONS']['positions']:
        element = position['element']
        coordinates = " ".join(position['coordinates'])
        atomic_positions.append(f"{coordinates}")
    #calculate total electron in this system
    pseudopotentials_directory = "/home/qianwang/workplace/pseudopotential/nc-sr-05_pbe_standard_upf"
    total_electrons = calculate_total_electrons(qe_input, pseudopotentials_directory)

    # #generate files which band calculation used.
    # generate_band_input(qe_input, total_electrons, args.num_kpoints)
    # generate_bands_input(prefix)
    # generate_projwfc_input(prefix)
    # generate_dos_input(prefix)

    # #run band calculation
    # print("Run scf calculation ...")
    # os.system(f"{args.pw_cmd} <scf.in > scf.out")
    # print("Done.\n")
    # print("Run dos calculation ...")
    # os.system('dos.x <dos.in > dos.out')
    # plot_dos(prefix)
    # print("Done.\n")
    # print("Run band calculation ...")
    # os.system(f"{args.pw_cmd} <band.in > band.out")
    # print("Done.\n")
    # print("Run bands calculation ...")
    # os.system("bands.x <bands.in > bands.out")
    # print("Done.\n")
    # print("Run projwfc calculation ...")
    # os.system("projwfc.x <projwfc.fat.in > projwfc.fat.out")
    # print("Done.\n")

    #extract fermi energy from scf.out file
    system, efermi = ('metal',  13.5578)
    if not efermi:
        print(f"Can not find fermi energy in {prefix}.save/data-file-schema.xml file.")
        sys.exit()

    #extract band energy in band.out file
    #band_energy is a two dimensional array, means: band_energy[kpoint][band].
    band_energy = extract_band_energy(prefix)
    band_energy_max = max(max(row) for row in band_energy) -efermi
    band_energy_min = min(min(row) for row in band_energy) -efermi

    #extract KPATH.in label and kpoints number
    label, nkpoints = extract_label("KPATH.in",args.num_kpoints)

    #extract projection orbital
    #projband is the weight for each atoms' orbital
    #projband is a three dimensional array, means: projband[angular number][kpoint][band]
    #It should used with atom_symbol and orbital_symbol
    atom_index, atom_symbol, orbital_symbol, projband = extract_orbital("fatband.projwfc_up")

    #write band energy and projection band file
    with open('projband.dat','w') as file:
        for i in range(len(atom_symbol) ):
            file.write(atom_symbol[i])
            file.write(f"       {orbital_symbol[i]} \n")
            for j in range(len(projband[i])):
                for k in range(len(projband[i][j])):
                    file.write(f"{band_energy[j][k]}        {projband[i][j][k]}\n")

    #find projector based on projband and efermi, then find the energy window
    projector_atom_index = []
    projector_atom_symbol = []
    projector_orbital = []
    with open("projector.dat",'w') as file:
        file.write(f"{system} \n")
        flag = False
        for i in range(len(atom_symbol)):
            for j in range(len(projband[i])):
                for k in range(len(projband[i][j])):
                    temp_energy = band_energy[j][k]
                    temp_projband = projband[i][j][k]
                    #judge the contribution for atom ortital to near fermi energy bands
                    if abs(temp_energy - efermi) < 3 and temp_projband > 0.1:
                        projector_atom_index.append(atom_index[i])
                        projector_atom_symbol.append(atom_symbol[i])
                        projector_orbital.append(orbital_symbol[i])
                        file.write(f"{atom_index[i]} : {atom_symbol[i]} : f={','.join(atomic_positions[int(atom_index[i])-1].strip().split())} : {orbital_symbol[i][1].lower()}\n")
                        flag = True
                        break
                if flag:
                    flag = False
                    break
        
        #find the dis_win_min and dis_win_max, which is the lowerest and highest band energy for projector
        #Note: In EPW, it doesn't has dis_win_min parameter, but dis_win_mim = dis_froz_min in most time. 
        #      So we output the dis_froz_min
        #at meanwhile, get the projector's weigh sum in all band.
        #at meanwhile, find the lowerest band index
        projector_all = np.zeros((len(projband[0]),len(projband[0][0])))
        dis_win_min = 9999.9
        dis_win_max = -9999.9
        band_min_index = 9999
        for l in range(len(projector_atom_symbol)):
            for i in range(len(atom_symbol)):
                if atom_index[i] == projector_atom_index[l] and atom_symbol[i] == projector_atom_symbol[l] and orbital_symbol[i] == projector_orbital[l]:
                    for j in range(len(projband[i])):
                        for k in range(len(projband[i][j])):
                            projector_all[j][k] += projband[i][j][k]
                            if projband[i][j][k] > 0.1:
                                dis_win_min = min(dis_win_min, band_energy[j][k])
                                dis_win_max = max(dis_win_max, band_energy[j][k])
                                band_min_index = min(band_min_index, k)
        #check dose the band[band_min_index] and band[band_min_index -1] has gap
        has_gap = True
        for j in range(len(band_energy)):
            if band_energy[j][band_min_index] - band_energy[j][band_min_index - 1] < 0.1:
                has_gap = False
                break
        if has_gap:
            file.write(f"exclude_bands = 1:{band_min_index}\n")
        file.write(f"dis_win_max = {dis_win_max+0.05}\n")
        file.write(f"dis_win_min = {dis_win_min- 0.05}\n")
        file.write(f"dis_froz_min = {dis_win_min-0.05}\n")

        #find the dis_froz_max, which is the non-projector minimum energy
        dis_froz_max = dis_win_max
        for j in range(len(projector_all)):
            for k in range(len(projector_all[j])):
                if band_energy[j][k] <= dis_win_max and band_energy[j][k] >= efermi and projector_all[j][k] < 0.5:
                    dis_froz_max = min(dis_froz_max, band_energy[j][k])
        file.write(f"dis_froz_max = {dis_froz_max - 0.05}\n")

    #write project_all data
    with open('projband_all.dat','w') as file:
        for i in range(len(projector_all)):
            for j in range(len(projector_all[0])):
                file.write(f"{band_energy[i][j]}        {projector_all[i][j]}\n")

    #through orbital_symbol, generate orbital_index
    orbital_index = []
    for i in range(len(orbital_symbol)):
        if 'S' in orbital_symbol[i]:
            orbital_index.append(0)
        elif 'P' in orbital_symbol[i]:
            orbital_index.append(1)
        elif 'D' in orbital_symbol[i]:
            orbital_index.append(2)
        elif 'F' in orbital_symbol[i]:
            orbital_index.append(3)
        else:
            print(f"Not support angular number : {orbital_symbol[i]}")
            sys.exit()

    #write qe_pp.py file which used to run dfttoolbox python script
    with open('qe_pp.py','w') as file:
        file.write("from DFTtoolbox.qe import postproc\n")
        file.write("import os\n")
        file.write("run_task=[1,2,3,4]\n")
        file.write("wkdir=os.path.dirname(os.path.realpath(__file__))\n")
        file.write(f"Ef= {efermi}\n")
        file.write(f"kdiv={nkpoints}\n")
        file.write(f"klabel={label}\n")
        file.write(f"Ebound=[{band_energy_min-5},{band_energy_max+5}]\n")
        file.write("state_grp=[")
        for i in range(len(atom_index)):
            file.write(f"['{atom_index[i]}:{atom_index[i]}/{orbital_index[i]}/a/a'],")
        file.write(f"['{atom_index[0]}:{atom_index[-1]}/a/a/a']")
        file.write("]\n")
        file.write('''# Main ================================================================
pp=postproc(wkdir)
for task in run_task:
    if task==1: #'band_read':
        pp.band_read(Ef=Ef,bandfile='band.out')
    elif task==2: #'band_plot':
        pp.band_plot(kdiv=kdiv,klabel=klabel,Ebound=Ebound)
    elif task==3: #'fatband_read':
        pp.fatband_read(Ef=Ef,projout='projwfc.fat.out',projprefix='fatband')
    elif task==4: #'fatband_plot':
        pp.fatband_plot(state_grp=state_grp,kdiv=kdiv,klabel=klabel,Ebound=Ebound)
    elif task==5: #pdos_read:
        pp.pdos_read(Ef=Ef)
    elif task==6:
        pp.pdos_plot(state_grp=state_grp,Ebound=Ebound)''')
    
    #run qe_pp.py script, please ensure your environment has dfttoolbox. You can install through pip install dfttoolbox
    subprocess.run(['python', 'qe_pp.py'])

