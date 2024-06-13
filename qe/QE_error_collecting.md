# <div align="center"> <span style="color:red"> QE 报错集合 </span> </div>

## 电声耦合计算的报错集合 
### <span style="color:green"> 1. 声子计算加速自洽  </span> </div>
1. alpha_mix(1)=0.3 一般用0.3， 算的比较快，值越大算的越慢。**所以我采用的默认参数alpha_mix(1)=0.3, 如果你自己需要高精度，自己调整**。
2. ph.in中的tr2_ph=1.0d-14开关设置14就很好了，有时候甚至可以设置13.
    
### <span style="color:green"> 2. 声子计算自洽不收敛  </span> </div>
有时候会出现如下的报错：**在计算某一个q点的不可约表示电子自洽不收敛。

```shell

     Atomic displacements:
     There are  197 irreducible representations # 这里表示总共有197个不可约表示
      
     Representation     1      1 modes -  To be done
      
     Representation     2      1 modes -  To be done
     ....
     ....
     ....
     Representation   196      1 modes -  To be done
 
     Representation   197      1 modes -  To be done
 
 
     Representation #   1 mode #   1 # 计算第一个不可约表示的自洽
 
     Self-consistent Calculation
 
      iter #   1 total cpu time :  1751.5 secs   av.it.:   9.3
      thresh= 1.000E-02 alpha_mix =  0.300 |ddv_scf|^2 =  8.194E-06
 
      iter #   2 total cpu time :  1883.1 secs   av.it.:  26.1
      thresh= 2.863E-04 alpha_mix =  0.300 |ddv_scf|^2 =  1.522E-03

    
      iter # 149 total cpu time : 16231.8 secs   av.it.:  22.0 
      thresh= 1.619E-06 alpha_mix =  0.300 |ddv_scf|^2 =  1.334E-10
          
      iter # 150 total cpu time : 16330.3 secs   av.it.:  23.8 
      thresh= 1.155E-06 alpha_mix =  0.300 |ddv_scf|^2 =  1.702E-10
    
End of self-consistent calculation  # 计算第一个不可约表示的自洽失败。没有收敛
    
No convergence has been achieved 
```


1. 调小alpha_mix(1), 默认值是0.7. 例如调小至0.01。在每次更新迭代自洽势时混入的旧势的比例**
2. 调低tr2_ph，默认值时1e-12, 通常用1e-14。它控制自洽的收敛精度
3. nmix_ph默认值是4，适当提高到8~20. 在混合势中的自洽迭代次数。可以显著加快收敛速度，但代价是使用更多内存。
4. niter_ph 默认值时100, 可以适当提高到500
5. 增加scffit.in和scf.in中的nbnd. 关于如何判断你的体系由多少个价带多少个空带，可以查看scffit.out中 "occupation numbers" 关键字下面的内容.
```shell
grep "occupation numbers" scffit.out
```
1. 有时候可能不是在第一个不可约表示因为电子自洽停下来不收敛的，那么可以使用ph.in中的以下开关从停下来的那个不可约表示续算。这是一种更加精确的续算方法。
```shell
start_irr=210
last_irr=240
```

### <span style="color:green"> 3. 声子虚频  </span> </div>
4. 自洽计算时的degauss是一个很重要的参数，一般用0.02，但是如果出现一个很小的虚频的话，可以试一试0.05

### <span style="color:green"> 4. 声子不可约表示错误  </span> </div>
报错 from set_irr_sym_new : error。wrong representation。在计算声子那一步时，在某些不可约q点计算电声耦合时，总是出现这样的一个错误。

解决方法:
```shell
#1. 进入qe的安装目录的PHono的PH模块中：
cd PHonon/PH/

# 2. 打开随机矩阵这个文件
vi random_matrix.f90

# 3. 取消这行代码的注释
!! #defineUNIFORM DISTRIB

修改后正确的是：
#defineUNIFORM DISTRIB
```

### <span style="color:green"> 5. 报错 q-mesh breaks symmetry  </span> </div>
这是因为读入结构精度的问题, 一定要保证分q点目录自洽文件和不分q点的自洽文件中晶格参数和原子坐标一模一样。

目前qe_main.py读入结构的晶格矩阵和坐标的方式有多种，下面贴出源码一一解答：
```python
# 方式一: 如果`-i 输入结构文件`是relax.out， 那么就将relax.out中的坐标和晶格按照字符串读进来
if self.input_file_path.name == "relax.out" and self.input_file_path.exists():
    ...
    ...
# 方式二: 如果`-i 输入结构文件`是POSCAR, CONTCAR, *.vasp， 那么就将其中的坐标和晶格按照字符串读进来
elif self.input_file_path.name == "POSCAR" or self.input_file_path.name == "CONTCAR" or ".vasp" in self.input_file_path.name:
    ...
    ...
# 方式三: 如果`-i 输入结构文件`是除上述提到的文件名命名的结构文件，比如:scffit.out, scffit.in, scf.out, scf.in, ph.in, ph.out, 那么就用pymatgen的方式读入该输入文件
else:
    ...
    ...
```

### <span style="color:green"> 6. 报错：Attempting to use an MPI routine before initializing MPICH。  </span> </div>

当你在执行/work/home/may/software/qe-7.1/bin/lambda.x <lambda.in > lambda.out命令时，出现了上述错误，推测原因可能是你的某一个q点的电声耦合文件elph.inp_lambda.*的DOS和EF值和其它文件的不一样，所以lambda.in计算失败导致。  

**需要你仔细检查是不是该q点计算时scffit.in的参数与其它q点的不一致导致。因为有时候我们在算某一个q点的时候scffit.in中的nbnd设置的不够大，导致该q点的计算不收敛，为了加速收敛，手动修改了其nbnd，导致算出来的该q点的DOS和EF与其它q点的有微小的不一致。**

### <span style="color:green"> 7. 报错 from davcio : error; error reading file /.../tmp/_ph0/La1Ce1Y1Th1Be4H32.q_2/La1Ce1Y1Th1Be4H32.wfc27  </span> </div>
 
如果你先用112核并行计算scffit和scf和split_ph.in，然后你因为某种原因声子计算中断了, 你在续算声子的时候先用40核重新做scffit和scf和scffit和scf和split_ph.in(recover=.true.), 就会报这个错。 

原因: 当使用112核并行时, 就会在tmp/_ph0/La1Ce1Y1Th1Be4H32.q_2/*wfc*中产生112个波函数文件. 当使用40核并行时, 就会在tmp/_ph0/La1Ce1Y1Th1Be4H32.q_2/*wfc*中产生40个波函数文件. QE在读取波函数文件时，会把全部的波函数文件都读进来, 所以用40核并行时会读入40个核计算过的wfc和112核计算的wfc.

解决方法：把tmp目录下所有的wfc都删除了
```shell

find ./ -name "La1Ce1Y1Th1Be4H32.save" | xargs rm -fr
find ./ -type f \( -name "*wfc*" -o -name "*dwf*"  -o -name "*prd*" -o -name "*bar*"  -o -name "*recover*"   -o -name "*mixd*"  \) | xargs rm -rf

```

### <span style="color:green"> 8. 报错 电声耦合因为内存原因停止。此时动力学矩阵计算已经完成。如何完成后续计算。  </span> </div>

首先将split_ph.in中的trans=.true.,改为trans=.false., 删除recover=.true.,

以下文件非常重要, 我们以La1Ce1Th2Be4H32体系的第2个q点为例

如果采用 ldisp=.true. 和 start_q=2 和 last_q=2 的方式计算电声耦合
```shell
La1Ce1Th2Be4H32.dyn2
tmp/_ph0/La1Ce1Th2Be4H32.q_2/La1Ce1Th2Be4H32.dv1
tmp/_ph0/La1Ce1Th2Be4H32.q_2/La1Ce1Th2Be4H32.dv1_paw1
tmp/_ph0/La1Ce1Th2Be4H32.phsave/patterns.{1..27}.xml 
# 总共27个q点，所以有27个patterns文件, 并且其中的内容有 <QPOINT_NUMBER>2</QPOINT_NUMBER>, 这里的2代表 ldisp=.true. 模式下 第2个q点的编号
```

如果采用 ldisp=.false. 和 XQ1 XQ2 XQ3 的方式计算电声耦合
```shell
La1Ce1Th2Be4H32.dyn2
tmp/_ph0/La1Ce1Th2Be4H32.q_2/La1Ce1Th2Be4H32.dv1
tmp/_ph0/La1Ce1Th2Be4H32.q_2/La1Ce1Th2Be4H32.dv1_paw1
tmp/_ph0/La1Ce1Th2Be4H32.phsave/patterns.1.xml
 # 总共27个q点，但是只显示一个patterns文件，并且编号为1, 并且其中的内容有 <QPOINT_NUMBER>1</QPOINT_NUMBER>, 这里1代表 ldisp=.false. 模式下所有的q点编号都是1
```

最后重新计算scffit.in 和 scf.in 和 split_ph.in