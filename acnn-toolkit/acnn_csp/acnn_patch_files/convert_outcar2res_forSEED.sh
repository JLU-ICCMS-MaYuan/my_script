#!/bin/bash

for task_dir in */; do
    # 去掉末尾斜杠
    task_name="${task_dir%/}"

    # 调试：显示真实任务名（包括隐藏字符）
    printf "Processing directory: '%s' → task_name='%s'\n" "$task_dir" "$task_name"

    outcar_path="$task_dir/OUTCAR"
    res_file="$task_dir/${task_name}.res"

    if [[ -f "$outcar_path" ]]; then
        # 执行 outcar2seed（抑制其 stdout？或者确认它不该输出到 stdout）
        if ! outcar2seed "$outcar_path" "$res_file"; then
            echo "Error: outcar2seed failed for $task_name" >&2
            exit 1
        fi

        # 使用 | 作为 sed 分隔符，避免 / 冲突
        if ! sed -i "s|cabal-in-out|$task_name|g" "$res_file"; then
            echo "sed failed for $res_file" >&2
            exit 1
        fi

        echo "Processed: ${task_name}.res -> Replaced 'cabal-in-out' with '$task_name'"
    else
        echo "  OUTCAR not found in $task_dir"
    fi
done