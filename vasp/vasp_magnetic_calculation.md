# VASP 磁性材料计算

## 必须要做的设置
```shell
ISPIN    = 2
LASPH    = .TRUE.
LORBIT   = 11
```

##  <span style="color:yellow">  U\J的设置

以CeSc2H24为例说明, 只设置U值，此时可以用LDAUTYPE = 2
```shell
LDAU     = .TRUE.
LDAUTYPE = 2 # LDAUTYPE=2意味着只需要设置U值即可，此时的U值是U-J的最终结果
LDAUL    = 3 -1 -1
LDAUU    = 4.0 0.0 0.0
LMAXMIX  = 6
```

有些论文中不仅给出了U值，还给出了J值，此时可以用LDAUTYPE=1
```shell
LDAUTYPE = 1 # LDAUTYPE=1意味着只需要设置U值即可
LDAUL    = 3 -1 -1
LDAUU    = 4.5 0.0 0.0
LDAUJ    = 0.5 0.0 0.0
LMAXMIX  = 6
```

##  <span style="color:yellow">  磁矩的设置

### 共线
做共线计算（ISPIN = 2），不需要 用 vasp_ncl。恰恰相反，推荐使用 vasp_std
```shell
ISYM     = 2 # 或者不设置ISYM直接采用默认值即可
MAGMOM   = 1*3  2*0  24*0
LSORBIT  =.FALSE. # 共线可以开SOC，也可以不开SOC，但是结构优化一般不开启SOC
```

### 非共线
非共线计算必须用vasp_ncl
```shell
ISYM     = -1 # 非共线必须关闭对称性
MAGMOM   = 0 0 3    6*0     72*0
GGA_COMPAT = .FALSE. # 提高GGA计算非共线的精度
LSORBIT =.TRUE. # 非共线必须开SOC
```
