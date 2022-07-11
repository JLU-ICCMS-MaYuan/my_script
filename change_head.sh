#!/bin/bash
for fpy in *.py; do
	echo $fpy
	head -n 2 $fpy	
	sed -i '1 i #!/work/home/may/miniconda3/bin/python3' $fpy
	sed -i '/\/mayuan\/miniconda3\//d' $fpy
	sed -i 's/mayuan/may/g' $fpy
	sed -i 's/xieyu/lhy/g' $fpy
done

