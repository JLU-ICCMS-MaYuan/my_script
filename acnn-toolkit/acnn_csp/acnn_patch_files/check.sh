#!/bin/bash

# 一个简化的只读脚本，用于分析和显示配置，不会修改任何文件。

# --- 定义文件路径 ---
relax_script="RELAX/dyn_batch_relax_lmp"
dft_incar="DFT/dyn_vasp_in"
dft_sub="DFT/sub.sh"
pot_sub="POT/sub.sh"
auto_script="auto" # 新增 auto 脚本路径定义

echo "--- 配置分析报告 ---"
echo "当前目录: $(pwd)"

# --- 1. 分析 RELAX 配置 ---
echo ""
echo "----------------------------------------"
echo "1. RELAX 分析 ($relax_script)"
echo "----------------------------------------"
if [ -f "$relax_script" ]; then
    echo "  关键变量:"
    # 提取并打印变量，保留注释
    grep -E '^(prj|types|press|frame|group|warp|job_max|JOB_NAME)=' "$relax_script" | while read -r line; do
        var_name=$(echo "$line" | cut -d'=' -f1)
        value_part=$(echo "$line" | cut -d'=' -f2-)
        # 修正：只清理 value_part, 而不是整个 line, 避免变量名重复
        cleaned_value=$(echo "$value_part" | sed 's/"//g')
        printf "    - %-8s = %s\n" "$var_name" "$cleaned_value"
    done

    echo ""
    echo "  SLURM 配置 (来自 TASK 函数):"
    # 修正：直接从文件中 grep, 不再使用复杂的 awk 匹配
    grep '#SBATCH' "$relax_script" | grep -E 'job-name|partition|time|exclude|ntasks-per-node' | sed 's/^/    - /'

    echo ""
    echo "  环境激活 (来自 TASK 函数):"
    # 修正：直接从文件中 grep, 并排除顶部的 in.seed
    grep 'source ' "$relax_script" | grep -v 'in.seed' | sed 's/^/    - /'
else
    echo "  文件未找到。"
fi

# --- 2. 分析 DFT 配置 ---
echo ""
echo "----------------------------------------"
echo "2. DFT 分析"
echo "----------------------------------------"
if [ -d "DFT" ]; then
    # 分析 DFT/dyn_vasp_in
    echo "  INCAR 参数 ($dft_incar):"
    if [ -f "$dft_incar" ]; then
        grep -E '^\s*(ENCUT|KSPACING)\s*=' "$dft_incar" | sed 's/^\s*/    - /'
    else
        echo "    文件未找到。"
    fi

    # 检查 POTCAR 文件
    echo ""
    echo "  POTCAR 文件 (在 DFT/ 目录下):"
    if ls DFT/POTCAR-* > /dev/null 2>&1; then
        echo "    - 状态: 已找到。"
    else
        echo "    - 状态: 未找到。"
    fi
    
    # 分析 DFT/sub.sh
    echo ""
    echo "  提交脚本 ($dft_sub):"
    if [ -f "$dft_sub" ]; then
        echo "    SLURM 配置:"
        grep '#SBATCH' "$dft_sub" | grep -E 'job-name|ntasks-per-node|time|partition|exclude' | sed 's/^/      - /'
        
        echo ""
        echo "    执行命令:"
        grep '^cmd=' "$dft_sub" | sed 's/^/      - /'

        echo ""
        echo "    环境激活:"
        grep 'source ' "$dft_sub" | sed 's/^/      - /'
    else
        echo "    文件未找到。"
    fi
else
    echo "  DFT 目录未找到。"
fi

# --- 3. 分析 POT 配置 ---
echo ""
echo "----------------------------------------"
echo "3. POT 分析"
echo "----------------------------------------"
if [ -d "POT" ]; then
    # 分析 POT/sub.sh
    echo "  提交脚本 ($pot_sub):"
    if [ -f "$pot_sub" ]; then
        echo "    SLURM 配置:"
        grep '#SBATCH' "$pot_sub" | grep -E 'job-name|ntasks-per-node|time|partition|exclude' | sed 's/^/      - /'
        
        echo ""
        echo "    环境激活:"
        grep 'source ' "$pot_sub" | sed 's/^/      - /'
    else
        echo "    文件未找到。"
    fi
else
    echo "  POT 目录未找到。"
fi

# --- 4. 分析自动化流程 (auto) --- [新增功能]
echo ""
echo "----------------------------------------"
echo "4. 自动化流程分析 (auto)"
echo "----------------------------------------"
if [ -f "$auto_script" ]; then
    echo "  自动化脚本: $auto_script"
    echo ""
    
    # 检查 RELAX 步骤是否使用 dyn_batch_relax_lmp
    echo "  RELAX 步骤检查:"
    # 使用 if/then 结构进行更可靠的检查
    if grep "dyn_batch_relax_lmp" "$auto_script" ; then
        echo "    - 状态: 确认使用 'dyn_batch_relax_lmp' 进行结构优化。"
    else
        echo "    - 状态: 未在 RELAX 步骤中找到 'dyn_batch_relax_lmp' 的明确调用。"
    fi

    # 检查其他主要执行步骤
    echo ""
    echo "  主要执行步骤:"
    grep -E 'cd (DFT|XSF|PD|POT|RELAX)' "$auto_script" | sed 's/cd/    - 进入目录:/'

else
    echo "  自动化脚本 'auto' 未找到。"
fi


echo ""
echo "----------------------------------------"
echo "分析完成。"
echo "----------------------------------------"