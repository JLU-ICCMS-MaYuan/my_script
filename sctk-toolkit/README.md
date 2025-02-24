# <div align="center"> **sctk-toolkit**  </div>

网络帖子: https://blog.csdn.net/lielie12138/article/details/122471410
官方教程: https://mitsuaki1987.github.io/sctk/en/_build/html/index.html

# 1 计算自洽
```shell
qe_main.py -i relax.out -j slurm -p PRESS scf -m mode=scf ecutwfc=80 ecutrho=960 kpoints_sparse='8 8 8'  degauss=0.02 execmd='mpirun -np 8' npool=4 queue=lhy
```

# 2 计算声子但是不计算电声耦合
```shell
qe_main.py -i relax.out -j bash sctk -m mode=nosplit  dyn0_flag=False EPC_flag=False search_sym=.false. trans=.true.  execmd='mpirun -np 8'  qpoints='4 4 4' queue=lhy
```

# 3 计算电声耦合
```shell
# 先用 fermi_velocity.x 获得 elph_nbnd_min elph_nbnd_max 的数值, 即:
./fermi_velocity.x < scf.in > fermit_v.out

# 特别注意sctk修改过的电声耦合计算要用到search_sym=.false. trans=.false. elph_nbnd_min=2 elph_nbnd_max=6
qe_main.py -i relax.out -j bash sctk -m mode=nosplit  dyn0_flag=False EPC_flag=True search_sym=.false. trans=.false. elph_nbnd_min=2 elph_nbnd_max=6 execmd='mpirun -np 8' npool=8  qpoints='4 4 4'
```

# 总: 下面的4,5,6,7,8,9,10都可以用一个代码来解决:
```shell
# k_automatic=True,  nscf的k点通过qe自动产生
# k_automatic=False, nscf的k点坐标通过指定k点产生
# nbnd 一定要手动设置比默认大一点, nscf, twin里面用的上
# sctk_all 包含了: nscf, twin, kel, lambda_mu_k, scdft_tc, deltaf, qpdos
qe_main.py -i relax.out -j bash -l debug sctk -m mode=sctk_all nbnd=20  occupations=tetrahedra_opt  execmd='mpirun -np 8' npool=8  qpoints='4 4 4' k_automatic=True queue=lhy
```

# 分: 单独计算4,5,6,7,8,9,10
## 4 计算非自洽nscf
```shell
# 直接再当前目录计算即可
# charge_density_dat和data_file_schema_xml路径不需要指定,直接默认使用 self.charge_density_dat = "tmp/H3S1.save/charge-density.dat" 和 self.data_file_schema_xml = "tmp/H3S1.save/data-file-schema.xml"
qe_main.py -i relax.out -j bash scf -m mode=nscf ecutwfc=80 ecutrho=960 nbnd=20 kpoints_dense='16 16 16'  degauss=0.02 execmd='mpirun -np 8' npool=8   k_automatic=True   occupations=tetrahedra_opt  queue=lhy
```

## 5 计算 twin
```shell
qe_main.py -i relax.out -j bash sctk -m mode=twin ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```

## 6 计算 kel
```shell
qe_main.py -i relax.out -j bash sctk -m mode=kel  ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```

## 7 计算 lambda_mu_k
```shell
qe_main.py -i relax.out -j bash sctk -m mode=lambda_mu_k  ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```

## 8 计算 scdft_tc
```shell
qe_main.py -i relax.out -j bash sctk -m mode=scdft_tc  ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```

## 9 计算 deltaf
```shell
qe_main.py -i relax.out -j bash sctk -m mode=deltaf  ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```

## 10 计算 qpdos
```shell
qe_main.py -i relax.out -j bash sctk -m mode=qpdos  ecutwfc=80 ecutrho=960 nbnd=20 qpoints='4 4 4' degauss=0.02 execmd='mpirun -np 8' npool=8 k_automatic=False queue=lhy
```