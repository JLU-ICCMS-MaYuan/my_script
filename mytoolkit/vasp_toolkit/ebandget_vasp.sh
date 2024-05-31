#!/bin/bash

echo "NOTE --------------------"
echo "    This script have to run in eband directory !"
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
echo "    Use vaspkit get the Band-Structure"
echo -e "21\n211" | vaspkit > \dev\null

echo "NOTE --------------------"
echo "    Use vaspkit get the Projected-Orbits-Band-Structure"
echo -e "21\n213" | vaspkit > \dev\null


echo "NOTE --------------------"
echo "    Use vaspkit get the Projected-Elements-Band-Structure"
echo -e "21\n215" | vaspkit > \dev\null

echo "NOTE --------------------"
echo "    Now you have got PBAND_ELEMENTS.dat, PBAND_*.dat BAND.dat"
