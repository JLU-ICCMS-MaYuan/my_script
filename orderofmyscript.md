超胞法（有限位移方法）声子
```shell
vasp_main.py -i CONTCAR-rvf/200.0/CONTCAR -w ./ -j slurm -p 200 phono -m mode=disp encut=800 kpoints="20 20 20" queue=lhy supercell="2 2 2" ismear=1
```

超胞法（有限位移方法）后处理
```shell
vasp_main.py -i POSCAR-init -w ./ data -m mode=dispprog supercell="2 2 2" mp="8 8 8" spectrum=False dos=True
```

qe结构弛豫
```shell
qe_main.py -i POSCAR -w ./ -j slurm -p 200 relax -m mode=relax-vc ecutwfc=80 ecutrho=960 forc_conv_thr=1.0d-6 etot_conv_thr=1.0d-8 queue=lhy press_conv_thr=0.001 kpoints_dense="20 20 20"
```

qe sciffit 自洽计算



mytoolkit 命令

将 cif 转化为 vasp
```shell
tool_main.py -i Ba3Si23.cif -w ./ convert -m vasp
```

