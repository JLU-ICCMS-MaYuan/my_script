# 脚本使用说明以及一些便利操作

### **1. 当你想获得所有成功优化的结构的焓值以及整个结构的目录名称时，比如这样的形式：**
```shell
POSCAR-0 -0.46534553
POSCAR-1 -0.46534984
POSCAR-2 -0.46348250
POSCAR-3 -0.46534907
POSCAR-4 -2.00337553
POSCAR-5 -0.46390904
POSCAR-6 -0.46381637
POSCAR-7 -0.46535080
POSCAR-8 -0.46535194
POSCAR-9 -0.46389392
```
### 操作步骤：
1. 获得所有优化成功的目录的位置
```shell
checkvasp_vasplog_converge.py | grep ok| awk '{print $1}' > check.log
```

2. 提取名称和焓值
```shell
for i in `cat check.log`; do cd $i; filename=$(basename "$i"); E=$(get_dH_peratoms.py | awk '{print $1}'); echo $filename $E >> ../dH.dat; done
```

3. 按照焓值排序
```shell
sort -n -k2 dH.dat -o dH.dat
# -n  表示按照数值排序而不是按照字典顺序排序
# -k2 表示按照第二列进行排序。
# -o  表示将修改好的内容导入dH.dat文件中
```