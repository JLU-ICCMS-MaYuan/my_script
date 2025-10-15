#!/usr/bin/env python3
import os
import subprocess
import sys
import re

def process_tasks(pressure):
    """
    遍历当前目录下的所有一级子目录，处理 OUTCAR 文件，
    并将生成的 .res 文件中的 'cabal-in-out' 替换为文件名。

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
                        check=True,
                        capture_output=True,
                        text=True
                    )
                    
                    # --- 整合自 Bash 脚本的替换功能 ---
                    
                    # 提取不带 .res 后缀的文件名，作为替换的新字符串
                    new_name = os.path.splitext(output_name)[0]
                    
                    try:
                        # 读取文件内容
                        with open(output_path, 'r', encoding='utf-8') as f:
                            content = f.read()

                        # 使用 re.sub() 进行替换
                        modified_content = re.sub(r'cabal-in-out', new_name, content)
                        
                        # 将修改后的内容写回文件
                        with open(output_path, 'w', encoding='utf-8') as f:
                            f.write(modified_content)
                        
                        print(f"Processed: {output_name} -> Replaced 'cabal-in-out' with '{new_name}'")

                    except FileNotFoundError:
                        print(f"  Error: Output file {output_path} not found.", file=sys.stderr)
                        sys.exit(1)
                    except Exception as e:
                        print(f"  Error processing file {output_path}: {e}", file=sys.stderr)
                        sys.exit(1)

                except subprocess.CalledProcessError as e:
                    print(f"  Error executing outcar2seed: {e.stderr}", file=sys.stderr)
                    sys.exit(1)
            else:
                print(f"  OUTCAR not found in {target_dir}")
        else:
            print(f"  Directory {target_dir} not found")

if __name__ == "__main__":
    try:
        if len(sys.argv) > 1:
            pressure = sys.argv[1]
        else:
            print("Usage: python script_name.py <pressure>")
            sys.exit(1)
    except (ValueError, IndexError):
        print("Invalid pressure value. Please provide an integer.")
        sys.exit(1)

    process_tasks(pressure)
