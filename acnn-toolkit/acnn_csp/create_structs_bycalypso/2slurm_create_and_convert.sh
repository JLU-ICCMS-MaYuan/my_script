#!/bin/bash
#SBATCH --job-name=caly_process_unfinished
#SBATCH --partition=multi-cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=1

# --- 可配置参数 ---
# 在这里设置您希望的最大并行任务数。
MAX_PARALLEL_JOBS=40
INPUT_FILE="unfinished.dat"

# --- 脚本主逻辑 ---

# 检查 unfinished.dat 文件是否存在且不为空
if [ ! -s "$INPUT_FILE" ]; then
    echo "错误: 未找到 $INPUT_FILE 文件，或文件为空。无法执行任务，正在退出。"
    exit 1
fi

echo "Starting job to process directories listed in $INPUT_FILE."
echo "Maximum parallel tasks: $MAX_PARALLEL_JOBS"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_NODELIST"

# 定义一个函数来执行每个子任务。
# 这个函数现在接收一个目录名作为参数。
run_task() {
    local dir_name=$1

    # 检查目录是否存在，以防万一
    if [ ! -d "$dir_name" ]; then
        echo "警告: 在 $INPUT_FILE 中列出的目录 '${dir_name}' 不存在，已跳过。"
        return
    fi

    echo "--> Starting task for ${dir_name}"
    # 目录已经由 Python 脚本创建，所以不再需要 'mkdir'
    cp create_and_convert.sh "${dir_name}"
    cd "${dir_name}"
    
    # 执行核心任务脚本
    sh create_and_convert.sh
    
    cd ../
    echo "<-- Finished task for ${dir_name}"
}


# --- ★★★ 核心修改：新的循环和并行控制 ★★★ ---
# 使用 while read 循环，逐行读取 unfinished.dat 文件中的目录名
while IFS= read -r dir_to_process || [ -n "$dir_to_process" ]; do
    # 1. 检查当前后台运行的任务数量
    while [ $(jobs -p | wc -l) -ge $MAX_PARALLEL_JOBS ]; do
        # 2. 如果任务数量达到上限，等待任意一个后台任务结束
        # 注意：您的原始脚本中使用 sleep 10。如果您希望更快地启动新任务，
        # 可以使用 wait -n (如果您的系统支持) 或者缩短 sleep 的时间。
        sleep 10
    done
    
    # 3. 启动新任务到后台，将从文件中读取的目录名作为参数传递给函数
    run_task "$dir_to_process" &

# 将 unfinished.dat 文件重定向为 while 循环的输入
done < "$INPUT_FILE"

# --- 清理阶段 ---
# 最后的 `wait` 命令至关重要，等待所有后台任务全部完成。
echo "All tasks have been launched. Waiting for the remaining tasks to complete..."
wait

echo "All tasks completed successfully at $(date)."
