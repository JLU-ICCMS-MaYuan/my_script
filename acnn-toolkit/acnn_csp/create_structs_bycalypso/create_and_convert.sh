#!/bin/bash

# --- Configuration ---
# 在这里设置你的体系名称和对应的元素
# 只需要修改这里，脚本的其他部分会自动更新
SYSTEM_NAME="LiPbH"
ELEMENTS="Li Pb H"

# --- Initial Cleanup ---
# 删除之前的结果和日志文件
# 注意：此处的 'rm -f' 会忽略不存在的文件，不会报错
rm -f -r results/
rm -f caly.log *py step

# --- Run Calculation ---
# 运行 calypso.x 并将输出重定向到日志文件
calypso.x > caly.log 2>&1

# --- Get Current Directory Name ---
# 用于最终文件的命名
a=$(basename "$(pwd)")

# 'wait' 命令在这里不是必需的，因为没有后台任务，但保留也无妨
wait

# --- Main Processing Loop ---
# 自动查找所有 POSCAR_* 文件并进行处理
echo "Starting processing for system: ${SYSTEM_NAME}"
for file in POSCAR_*; do
    # 如果找不到任何匹配的文件，则跳过循环以避免错误
    [ -e "$file" ] || continue

    # 从文件名中提取数字 (例如, 从 "POSCAR_123" 得到 "123")
    i=${file#POSCAR_}

    echo "Processing file: $file (Index: $i)"

    # 1. 将元素行插入到临时的 POSCAR 文件中
    #    注意这里使用了双引号，以便 ${ELEMENTS} 变量能被正确识别
    sed "6 i ${ELEMENTS}" "$file" > "POSCAR-$i"

    # 2. 运行第一个 cabal 命令，使用体系名变量
    cabal poscar res < "POSCAR-$i" > "../../Base/${SYSTEM_NAME}-${i}-${a}-N64.res"

    # 3. 运行 rres 命令
    rres "../../Base/${SYSTEM_NAME}-${i}-${a}-N64.res" -s "${SYSTEM_NAME}-${i}-${a}-N64-st"

    # 4. 运行第二个 cabal 命令
    cabal res res 0 < "../../Base/${SYSTEM_NAME}-${i}-${a}-N64.res" > "../../Base/${SYSTEM_NAME}-${i}-${a}-N64-st.res"

    # 5. 清理中间文件
    rm "../../Base/${SYSTEM_NAME}-${i}-${a}-N64.res"
    rm "POSCAR-$i"

done

echo "Script finished."