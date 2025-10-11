#!/bin/bash

# 一个简化的只读脚本，用于分析和显示配置，不会修改任何文件。

# --- 定义文件路径 ---
relax_script="RELAX/dyn_batch_relax_lmp"
dft_incar="DFT/dyn_vasp_in"
dft_sub="DFT/sub.sh"
pot_sub="POT/sub.sh"
auto_script="auto"

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
        cleaned_value=$(echo "$value_part" | sed 's/"//g')
        printf "    - %-8s = %s\n" "$var_name" "$cleaned_value"
    done

    echo ""
    echo "  SLURM 配置 (来自 TASK 函数):"
    grep '#SBATCH' "$relax_script" | grep -E 'job-name|partition|time|exclude|ntasks-per-node' | sed 's/^/    - /'

    echo ""
    echo "  环境激活 (来自 TASK 函数):"
    grep 'source ' "$relax_script" | grep -v 'in.seed' | sed 's/^/    - /'

    # --- [新增功能] 平台兼容性检查 ---
    echo ""
    echo "  平台兼容性提示:"
    # 检查是否包含特定于 Inspur 服务器的 LD_LIBRARY_PATH 设置
    if grep -q 'export LD_LIBRARY_PATH=/work/home/mayuan/bin:$LD_LIBRARY_PATH' "$relax_script"; then
        echo "    - Inspur 平台: 已找到推荐的 LD_LIBRARY_PATH 设置。"
    else
        echo "    - Inspur 平台: [注意] 如果在 Inspur 服务器上运行，建议在激活环境(source)附近添加以下行："
        echo "                   export LD_LIBRARY_PATH=/work/home/mayuan/bin:\$LD_LIBRARY_PATH"
    fi
    # --- [新增功能结束] ---

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

# --- 4. 分析自动化流程 (auto) ---
echo ""
echo "----------------------------------------"
echo "4. 自动化流程分析 (auto)"
echo "----------------------------------------"
if [ -f "$auto_script" ]; then
    echo "  自动化脚本: $auto_script"
    echo ""
    
    # 检查 RELAX 步骤是否使用 dyn_batch_relax_lmp
    echo "  RELAX 步骤检查:"
    if grep -q "cd RELAX" "$auto_script" && grep "dyn_batch_relax_lmp" "$auto_script" | grep -q "cd RELAX" -B 1; then
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