#!/usr/bin/env python3
import os
import subprocess
import sys

def process_tasks(pressure):
    """
    遍历当前目录下的所有一级子目录，并处理 OUTCAR 文件。

    Args:
        pressure (int): 计算所用的压力值。
    """
    current_dir = os.getcwd()
    
    # 遍历当前目录下的所有一级子目录
    for task_name in os.listdir(current_dir):
        task_path = os.path.join(current_dir, task_name)
        
        # 只处理目录
        if not os.path.isdir(task_path):
            continue


        target_dir = os.path.join(task_path, str(pressure))
        
        # 检查目标目录是否存在
        if os.path.isdir(target_dir):
            outcar_path = os.path.join(target_dir, "OUTCAR")
            
            # 检查 OUTCAR 文件是否存在
            if os.path.isfile(outcar_path):
                # 构造输出文件名
                output_name = f"{task_name}_{pressure}GPa.res"
                output_path = os.path.join(task_path, output_name)
                # 执行 outcar2seed 脚本
                try:
                    subprocess.run(
                        ["outcar2seed", outcar_path, output_path],
                        check=True,  # 如果命令失败，会抛出异常
                        capture_output=True,
                        text=True
                    )
                    print(f"Processing: {task_name}. Found OUTCAR. Converting to {output_name}...")
                except subprocess.CalledProcessError as e:
                    print(f"  Error executing outcar2seed: {e.stderr}", file=sys.stderr)
                    sys.exit(1)
            else:
                print(f"  OUTCAR not found in {target_dir}")
        else:
            print(f"  Directory {target_dir} not found")

if __name__ == "__main__":
    # 可以将压力值作为命令行参数传递
    # 例如：python get_res.py 50
    pressure = sys.argv[1]
    process_tasks(pressure)
