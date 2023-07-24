#!/bin/bash

echo "NOTE --------------------"
echo "    This script have to run in eledos directory !"
# 设置新的费米能级
new_fermi_energy=`grep E-fermi ../scf/OUTCAR | awk '{print $3}'`
echo "efermi =" $new_fermi_energy

# 打开DOSCAR文件
doscar_file="DOSCAR"
doscar_data=$(cat $doscar_file)

# 获取旧费米能级
old_fermi_energy=$(echo "$doscar_data" | head -n 6 | tail -n 1 | awk '{print $4}')

# 替换费米能级
new_doscar_data=$(echo "$doscar_data" | sed "s/$old_fermi_energy/$new_fermi_energy/g")

# 将修改后的数据写回DOSCAR文件
echo "$new_doscar_data" > $doscar_file

# 将替换好的DOCAR的前6行显示出来
echo "NOTE --------------------"
echo "    The final result show you !!"
echo "$new_doscar_data" | head -n 6 | tail -n 1 

echo "NOTE --------------------"
echo "    Use vaspkit get the PDOS for every element and TDOS"
echo -e "11\n113" | vaspkit > \dev\null
echo -e "11\n111" | vaspkit > \dev\null

echo "NOTE --------------------"
echo "    Get TDOS(Ef). The unit is states/eV/f.u. in VASP "
echo "                  (The unit is states/spin/Ry/f.u. in Quantum Espresso)"
echo "                  1 Ry = 13.605693122994 eV"
echo "                  1 eV = 0.0734986443513 Ry"
echo "                  DOS_vasp_qe = 0.0734986443513 * DOS_vasp_value"
echo "                  "
