#!/bin/bash

# 定义系统名称变量, 你需要修改的唯一参数
system_name="H12La2B1"

# 定义要检查的文件
files=("scffit.in" "scf.in" "s5_PhAssignQ.sh" "split_ph.in")

# 遍历每个文件并检查是否存在
for file in "${files[@]}"; do
    if [[ -f $file ]]; then
        echo "$file exists."
    else
        echo "$file does not exist."
    fi
done

# 生成 qlist.dat 文件
sed '1,2d' ${system_name}.dyn0 > qlist.dat
nq=`sed -n 2p ${system_name}.dyn0`

# 创建目录并复制文件
for j in `seq 1 $nq`; do 
    mkdir $j
    cp inter_dyn_${j} ${j}/${system_name}.dyn
done

# 处理每个 Q 点
for Q in `seq 1 $nq`; do
    cp scffit.in scf.in s5_PhAssignQ.sh $Q
    Q1=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $1 }')
    Q2=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $2 }')
    Q3=$(head ./qlist.dat -n${Q} | tail -n1 | awk '{ print $3 }')
    cat ./split_ph.in | sed -e "s/XQ1/$Q1/g" \
    | sed -e "s/XQ2/$Q2/g" \
    | sed -e "s/XQ3/$Q3/g" \
    > $Q/split_ph.in
done

for i in `seq 1 $nq`; do
    # 创建目录结构
    mkdir -p "${i}/tmp/_ph0/${system_name}.phsave"
    
    # ========== 修改后的复制dvscf文件部分 ==========
    # 检查源文件是否存在
    if ls "../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.dvscf"* 1>/dev/null 2>&1; then
        cp -v "../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.dvscf"* "${i}/tmp/_ph0/"
    else
        echo "错误: ../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.dvscf* 文件不存在"
    fi
    
    # ========== 复制patterns文件部分 ==========
    if [ -f "../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.phsave/patterns.1.xml" ]; then
        cp -v "../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.phsave/patterns.1.xml" \
              "${i}/tmp/_ph0/${system_name}.phsave/"
    else
        echo "错误: ../2.interpolation/2.fine/${i}/tmp/_ph0/${system_name}.phsave/patterns.1.xml 文件不存在"
    fi
    
    echo "已处理文件夹 ${i}"
    echo "----------------------------------------"
done

