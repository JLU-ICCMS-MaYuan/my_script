## 1. 如何设置原子间距离

1. 替换已知结构并优化，测量键长*0.7
2. POTCAR中的RCORE\*0.529177\*0.7
3. 设置1.0

## 2. 如何设置晶胞体积

1. 定组分结构预测，体积为晶胞的体积，晶胞体积的计算需要结合RCORE的半径计算体积并乘1.3
2. 设置为0
3. 变组分结构预测的情况下，体积为每原子的体积。

## 3. calypso变组分结构预测报错
```shell
  _______  _______  _              _______  _______  _______ 
 (  ____ \(  ___  )( \   |\     /|(  ____ )(  ____ \(  ___  )
 | (    \/| (   ) || (   ( \   / )| (    )|| (    \/| (   ) |/20
 | |      | (___) || |    \ (_) / | (____)|| (_____ | |   | |
 | |      |  ___  || |     \   /  |  _____)(_____  )| |   | |
 | |      | (   ) || |      ) (   | (            ) || |   | |
 | (____/\| )   ( || (____/\| |   | )      /\____) || (___) |
 (_______/|/     \|(_______/\_/   |/       \_______)(_______)
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
mkdir: cannot create directory ‘step_1’: File exists
cp: cannot stat ‘1’: No such file or directory
cp: cannot stat ‘2’: No such file or directory
cp: cannot stat ‘3’: No such file or directory
cp: cannot stat ‘4’: No such file or directory
cp: cannot stat ‘5’: No such file or directory
forrtl: severe (408): fort: (2): Subscript #1 of the array ATOMICTYPE has value 9 which is greater than the upper bound of 8

Image              PC                Routine            Line        Source             
calypso.x          00000000016FEAFF  Unknown               Unknown  Unknown
calypso.x          000000000070789B  readvasp_                  86  ReadVasp.F90
calypso.x          00000000006E133C  readoptdata_               41  ReadOptData.F90
calypso.x          00000000009000EE  MAIN__                     96  Main.F90
calypso.x          0000000000404692  Unknown               Unknown  Unknown
calypso.x          00000000017C9E20  Unknown               Unknown  Unknown
calypso.x          000000000040457A  Unknown               Unknown  Unknown

```

解决方案：

```shell
NumberOfFormula = 1 4 改成 NumberOfFormula = 1 1

如果你想做变胞的结构预测，比如做2倍胞的结构预测就用：NumberOfFormula = 2 2 
```


## 4. 定组分结构预测报错

```shell
(cage) [may@ln01 init]$ /work/home/may/workplace/5.calypso/35.Ce-Sc-H/3.calypso-acnn/calypso.x.acnn_customized 
forrtl: severe (24): end-of-file during read, unit -5, file Internal List-Directed Read
Image              PC                Routine            Line        Source             
calypso.x.acnn_cu  0000000001540ABA  Unknown               Unknown  Unknown
calypso.x.acnn_cu  00000000015400E4  Unknown               Unknown  Unknown
calypso.x.acnn_cu  00000000004124C2  readfile_                 776  ReadFile.F90
calypso.x.acnn_cu  00000000008DF18A  inirun_                    27  IniRun.F90
calypso.x.acnn_cu  00000000008E257C  MAIN__                     45  Main.F90
calypso.x.acnn_cu  0000000000402ABD  Unknown               Unknown  Unknown
calypso.x.acnn_cu  00000000015BE8A0  Unknown               Unknown  Unknown
calypso.x.acnn_cu  000000000040299E  Unknown               Unknown  Unknown

```


```shell
# 这个开关非常坑，如果做变组分预测，没有设置得话会报错，
如果做定组分预测设置了得话，会报错。所以做定组分结构预测一定删掉这一行
VSCEnergy= 0 0 
```

