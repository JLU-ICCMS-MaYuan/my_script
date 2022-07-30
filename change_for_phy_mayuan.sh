#!/bin/bash
for fpy in *.py; do
	echo $fpy
	echo "before change"
	head -n 2 $fpy	
	sed -i '/\/may\/miniconda3\//d' $fpy
	sed -i '1 i #!/public/home/mayuan/miniconda3/envs/cage/bin/python3' $fpy
	echo "after change"
	head -n 2 $fpy	
	sed -i 's/xieyu/lhy/g' $fpy

done
#qeSuperconductTc.py 119
#PP = os.path.abspath("/public/home/mayuan/POT/qe-pp/all_pbe_UPF_v1.5")

