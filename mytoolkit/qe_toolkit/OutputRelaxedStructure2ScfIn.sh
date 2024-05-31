#!/bin/bash


# 获得两个行号
CELL_PARAMETERS=$(grep -n 'CELL_PARAMETERS' scf.in | awk -F : '{print $1}')
K_POINTS=$(grep -n 'K_POINTS' scf.in | awk -F : '{print $1}')

# 获得这两个行好所在行的前一行
insert_position=$(expr $CELL_PARAMETERS - 1)
stop_delete_position=$(expr $K_POINTS - 1)

# 获得总原子数
Nat=$(grep 'number of k points' -B 2 relax.out |head -n 1|awk {'print($1)'})
# 获得晶格矩阵参数
StruLine=$(expr $Nat + 5)
grep 'CELL_' -A $StruLine relax.out |tail -n "$(expr $StruLine + 1)" > new_structure.out

#原地置换，原文件删除空行
sed  -i '/^$/d' new_structure.out 
# 其它模板
# grep -v '^$' filename #打印非空行
# sed  '/^$/d' filename #打印非空行
# sed  -i '/^$/d' filename #原地置换，原文件删除空行
# awk '!/^$/{print}' filename #打印非空行

# 注意这里不能直接使用 echo ${structinfo}， 这样会导致echo出来的值没有回车 shell 输出其他命令的结果时，需要保持原格式，比如不需要去除掉换行，加个 "" 包裹输出对象即可
# echo -e "$structinfo"

# 单引号会消除$取值符号的特殊含义，$只会被解析为字符本身
sed -i "${CELL_PARAMETERS}, ${stop_delete_position}d" scf.in
# sed "${K_POINTS}i ${structinfo}" scf.in
sed -i "${insert_position}r new_structure.out" scf.in
sed -i "${insert_position}r new_structure.out" scffit.in
