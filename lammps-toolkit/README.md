# LAMMPS 使用注意事项

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

data