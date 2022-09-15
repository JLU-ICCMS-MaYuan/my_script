超胞法（有限位移方法）声子
```shell
vasp_main.py -i CONTCAR-rvf/200.0/CONTCAR -w ./ -j slurm -p 200 phono -m mode=disp encut=800 kpoints=[20,20,20] queue=lhy supercell=[2,2,2] ismear=1
```

超胞法（有限位移方法）后处理
```shell
vasp_main.py -i ./POSCAR-init -w ./ data -m mode=dispprog supercell=[2,2,2]
```

mytoolkit 命令

将 cif 转化为 vasp
```shell
tool_main.py -i Ba3Si23.cif -w ./ convert -m vasp
```

