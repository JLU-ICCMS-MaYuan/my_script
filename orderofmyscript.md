# VASP 计算
超胞法（有限位移方法）声子
```shell
vasp_main.py -i CONTCAR-rvf/200.0/CONTCAR -w ./ -j slurm -p 200 phono -m mode=disp encut=800 kpoints="20 20 20" queue=lhy supercell="2 2 2" ismear=1 ncore=4
```

超胞法（有限位移方法）后处理
```shell
vasp_main.py -i POSCAR-init -w ./ data -m mode=dispprog supercell="2 2 2" mp="8 8 8" spectrum=True dos=True
```

DFPT（密度泛函微软理论）声子
```shell
vasp_main.py -i CONTCAR-rvf/200.0/CONTCAR -w ./ -j slurm -p 200 phono -m mode=dfpt encut=800 kpoints="20 20 20" queue=lhy supercell="2 2 2" ismear=1
```
DFPT（密度泛函微扰理论）后处理
```shell
vasp_main.py -i POSCAR-init -w ./ data -m mode=dfptprog supercell="2 2 2" mp="8 8 8" spectrum=True dos=True
```


# QE 计算

准备输入文件
```shell
relax -m mode=relax-vc core=64 npool=8
```

进行批量计算
```shell
prepare -m mode="relax-vc scffit scf nosplit" dyn0_flag=True core=28 npool=4 queue=liuhy &
prepare -m mode="relax-vc scffit scf nosplit" dyn0_flag=True core=64 npool=8 queue=normal &
```


qe结构弛豫在上一步结构弛豫的基础上做新的结构弛豫
```shell
qe_main.py -i ./POSCAR          -w ./      -p 200 -j bash relax -m mode=relax-vc core=4 npool=1 queue=local
```
```shell
qe_main.py -i relax.out -w ./ -j slurm -p 200 relax -m mode=relax-vc ecutwfc=80 ecutrho=960 forc_conv_thr=1.0d-6 etot_conv_thr=1.0d-8 queue=lhy press_conv_thr=0.001 kpoints_dense="20 20 20"
```

qe-scffit 自洽计算(稠密)
```shell
qe_main.py -i relax.out -w ./ -j slurm -p 200 scf -m mode=scffit queue=lhy forc_conv_thr=1.0d-5 etot_conv_thr=1.0d-7 ecutwfc=80 ecutrho=960 smearing=methfessel-paxton conv_thr=1.0d-9 mixing_beta=0.8 kpoints_dense='24 24 24'
```


qe scf 自洽计算(稀疏)
```shell
qe_main.py -i relax.out -w ./ -j slurm -p 200 scf -m mode=scf queue=lhy forc_conv_thr=1.0d-5 etot_conv_thr=1.0d-7 ecutwfc=80 ecutrho=960 smearing=methfessel-paxton conv_thr=1.0d-9 mixing_beta=0.8 kpoints_sparse='12 12 12'
```

qe phono (不分q点计算声子)
```shell
qe_main.py -i relax.out -w ./ -j slurm -p 200 phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=False queue=lhy
```

qe phono (q点计算声子)
```shell
qe_main.py -i relax.out -w ./ -j slurm -p 200 phono -m mode=nosplit qpoints='6 6 6' dyn0_flag=True queue=xieyu
qe_main.py -i ./relax.out -w ./ -p 压力值  -j pbs phono -m mode=split_from_dyn0 qpoints='6 6 6' 
```

qe 合并声子文件
```shell
qe_main.py -i ./relax.out -w ./ -p 300  -j pbs phono -m mode=merge
```



mytoolkit 命令

将 cif 转化为 vasp
```shell
tool_main.py -i Ba3Si23.cif -w ./ convert -m vasp
```




    
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash scf   -m mode=scffit core=4 npool=1 queue=local
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash scf   -m mode=scf core=4 npool=1 queue=local
    # 非自洽计算: 用于计算dos, 使用四面体方法, 可以获得没有展宽的dos结果, 更加准确
    qe_main.py -i ./relax.out -w ./test  -p 压力值 -j pbs scf   -m mode=nscf  core=1 npool=1 queue=local kpoints_dense="20 20 20"

    # no split method
    # 进行不分q点计算
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash phono -m mode=nosplit core=1 npool=1 queue=local qpoints='4 4 4' dyn0_flag=False 
        # default value : dyn0_flag=False
    # 只计算`*dyn0`
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash phono -m mode=nosplit core=4 npool=1 queue=local qpoints='4 4 4' dyn0_flag=True
    # spilit method 1
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=split_dyn0 core=1 npool=1 queue=local
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=merge core=1 npool=1 queue=local
    # split method 2
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash phono -m mode=split_assignQ core=1 npool=1 queue=local

    # 声子计算其它可使用的参数
    el_ph_nsigma=50  

    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  phono -m mode=q2r        core=1 npool=1 queue=local
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  phono -m mode=matdyn     core=1 npool=1 queue=local qinserted=50
    qe_main.py -i ./200.0/relax.out -w ./200.0/ -p 200 -j bash  dos   -m mode=dos core=1 npool=1 queue=local qpoints='8 8 8' ndos=500 

    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash   sc    -m mode=McAD       core=1 npool=1 queue=local top_freq=80 deguass=0.5 screen_constant=0.1 smearing_method=1
    qe_main.py -i ./200.0/relax.out -w ./200.0 -p 200 -j bash   sc    -m mode=eliashberg core=1 npool=1 queue=local temperature_points=10000 a2F_dos=a2F.dos3