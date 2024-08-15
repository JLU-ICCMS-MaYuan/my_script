# 流程大纲

帖子
```shell
https://tiandijunhao.github.io/2020/06/28/tdep-ji-suan-sheng-zi-pu/

https://tiandijunhao.github.io/2020/06/24/tdep-an-zhuang/

https://ollehellman.github.io/page/workflows/minimal_example_2.html

https://tiandijunhao.github.io/2020/06/24/blas-lapack-an-zhuang/

https://tdep-developers.github.io/tdep/

https://mixzeng.github.io/2020/12/27/tdep-install/

# 我的帖子
https://zhuanlan.zhihu.com/p/714272645
```

## <span style="color:red">  准备文件

```shell
phonopy --symmetry -c POSCAR
cp PPOSCAR POSCAR-unitcell # 或者 cp BPOSCAR POSCAR-unitcell
cp POSCAR-unitcell infile.ucposcar
generate_structure -d value#1 vaule#2 vaule#3
cp outfile.ssposcar POSCAR
cp outfile.ssposcar infile.ssposcar
```

## <span style="color:red"> 提取OUTCAR

```shell
python process_outcar_5.3.py OUTCAR --skip=100000
```

## <span style="color:red"> 提取力常数

提取获得力常数后，会产生一个输出文件`outfile.forceconstant`

```shell
extract_forceconstants -rc2 5 
```

## <span style="color:red"> 创建力常数文件符号链接`infile.forceconstant`

```shell
ln -s outfile.forceconstant infile.forceconstant
```

## <span style="color:red"> 准备声子谱输入文件`infile.qpoints_dispersion`

```shell
CUSTOM
50
6
   0.0000   0.0000   0.0000    0.5000   0.0000   0.5000     GM    X
   0.5000   0.0000   0.5000    0.6250   0.2500   0.6250     X     U
   0.3750   0.3750   0.7500    0.0000   0.0000   0.0000     K     GM
   0.0000   0.0000   0.0000    0.5000   0.5000   0.5000     GM    L
   0.5000   0.5000   0.5000    0.5000   0.2500   0.7500     L     W
   0.5000   0.2500   0.7500    0.5000   0.0000   0.5000     W     X  
```

## <span style="color:red"> 绘制声子谱 

导出数据
```shell
phonon_dispersion_relations -rp 
#产生输出文件：
#   outfile.dispersion_relations.gnuplot
#   outfile.dispersion_relations.hdf5
#   outfile.dispersion_relations、
#   outfile.group_velocities.gnuplot
#   outfile.group_velocities
#   outfile.mode_activity.csv
#   outfile.irrifc_secondorder
```

修改outfile.dispersion_relations.gnuplot

```shell
 set terminal png   size 500,350 enhanced font "CMU Serif,10"
 set output "dispersion.png"
```

出图
```shell
gnuplot --persist outfile.dispersion_relations.gnuplot
```