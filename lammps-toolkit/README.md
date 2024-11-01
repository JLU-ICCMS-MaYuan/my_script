# LAMMPS 使用注意事项

网络热帖：
1. https://wap.sciencenet.cn/blog-3437453-1264601.html?mobile=1
2. https://casea1998.cn/2022/12/09/Lammps1/
3. https://lammpscnv2.vercel.app/
4. 


input通常分为四部分：

```shell
Initialization
    units
    dimension
    boundary
    atom_style
    bond_style
    pair_style
    dihedral_style

System definition
    read_data
    read_restart
    lattice, region, create_box, create_atoms

Settings
    pair/angle/dihedral/improper_coeff
    special_bonds
    neighbor
    fix
    compute

Running
    minimize
    run
```

运行lammps
```shell
lmp -in input.lammps
```

编译vasp couple lammps: https://zhuanlan.zhihu.com/p/4018186574