#!/bin/bash

Nat=$(grep 'number of k points' -B 2 relax.out |head -n 1|awk {'print($1)'})
StruLine=$(expr $Nat + 5)
grep 'CELL_' -A $StruLine relax.out |tail -n `expr $StruLine + 1` > new_structure.out
sed  -i '/^$/d' new_structure.out


cell_parameters=$(grep -n 'CELL_PARAMETERS' relax.in | awk -F : '{print $1}')
k_points=$(grep -n 'K_POINTS' relax.in | awk -F : '{print $1}')
insert_position=$(expr $cell_parameters - 1)
stop_delete_position=$(expr $k_points - 1)
sed -i "${cell_parameters}, ${stop_delete_position}d" relax.in
sed -i "${insert_position}r new_structure.out" relax.in

