#!/usr/bin/env python3

import os
import sys
import time
import argparse
from multiprocessing import Pool
from tqdm import tqdm

# --- 默认配置区 ---
# 这些值现在是 argparse 参数的默认值
DEFAULT_BASE_DIRECTORY = 'Base'
DEFAULT_OUTPUT_FILE = 'res_files_deleted_log.txt'
DEFAULT_NUM_CORES = 8
# --- 配置区结束 ---

def find_res_files(root_dir):
    """
    递归查找指定目录下所有以 .res 结尾的文件。
    
    :param root_dir: 要搜索的根目录
    :return: 包含所有 .res 文件绝对路径的列表
    """
    res_file_paths = []
    print(f"正在从 '{os.path.abspath(root_dir)}' 目录中查找所有 .res 文件...")
    # 为了能实时看到查找进度，这里也用一个简单的 walk 计数
    # 注意：这会稍微减慢查找速度，因为需要先遍历一遍获取总目录数
    try:
        total_dirs = sum(1 for _ in os.walk(root_dir))
        with tqdm(os.walk(root_dir), total=total_dirs, desc="[1/3] 正在扫描目录", unit="dir") as pbar:
            for dirpath, _, filenames in pbar:
                for filename in filenames:
                    if filename.lower().endswith('.res'):
                        full_path = os.path.join(dirpath, filename)
                        res_file_paths.append(os.path.abspath(full_path))
    except Exception as e:
        print(f"\n错误：在扫描目录时发生问题: {e}")
        sys.exit(1)
        
    return res_file_paths

def delete_file(filepath):
    """
    删除单个文件的函数，用于并行处理。
    """
    try:
        os.remove(filepath)
    except Exception:
        # 在并行处理中，忽略单个文件的错误以避免屏幕输出混乱
        pass

def main():
    """
    主执行函数
    """
    # 1. 设置和解析命令行参数
    parser = argparse.ArgumentParser(
        description="在 'Base' 目录下批量查找并删除 .res 后缀的文件。",
        formatter_class=argparse.RawTextHelpFormatter # 保持帮助信息格式
    )
    parser.add_argument(
        '-l', '--limit',
        type=int,
        default=None,
        help='指定最多删除多少个文件。\n如果省略此参数，将删除所有找到的文件。'
    )
    parser.add_argument(
        '-c', '--cores',
        type=int,
        default=DEFAULT_NUM_CORES,
        help=f'指定用于删除文件的并行进程数。\n默认值为: {DEFAULT_NUM_CORES}。'
    )
    args = parser.parse_args()

    # 2. 检查 Base 目录是否存在
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = os.path.join(script_dir, DEFAULT_BASE_DIRECTORY)

    if not os.path.isdir(target_dir):
        print(f"错误：在脚本所在位置未找到 '{DEFAULT_BASE_DIRECTORY}' 目录。")
        print(f"预期目录的绝对路径是: {target_dir}")
        sys.exit(1)

    # 3. 查找所有 .res 文件
    all_found_files = find_res_files(target_dir)
    total_found = len(all_found_files)

    if total_found == 0:
        print("\n扫描完成，未找到任何 .res 文件。程序退出。")
        sys.exit(0)

    print(f"\n扫描完成！共找到 {total_found} 个 .res 文件。")

    # 4. 根据 --limit 参数确定最终要处理的文件列表
    if args.limit is not None and args.limit >= 0:
        print(f"用户指定最多删除 {args.limit} 个文件。")
        # Python 的切片操作会自动处理 limit > total_found 的情况
        files_to_process = all_found_files[:args.limit]
    else:
        files_to_process = all_found_files
    
    num_to_delete = len(files_to_process)
    if num_to_delete < total_found:
        print(f"根据限制，将实际处理其中的 {num_to_delete} 个文件。")

    # 5. 将待删除文件路径写入日志文件
    try:
        with open(DEFAULT_OUTPUT_FILE, 'w', encoding='utf-8') as f:
            for path in tqdm(files_to_process, desc="[2/3] 正在将路径存入日志", unit="file"):
                f.write(path + '\n')
        print(f"所有待删除文件路径已记录到 '{DEFAULT_OUTPUT_FILE}' 文件中。")
    except IOError as e:
        print(f"错误：无法写入日志文件 '{DEFAULT_OUTPUT_FILE}'. 错误信息: {e}")
        sys.exit(1)

    # 6. 执行删除操作
    print("\n" + "!"*20)
    print(f"警告：将在 3 秒后使用 {args.cores} 个进程开始删除 {num_to_delete} 个文件...")
    print("按 Ctrl+C 可紧急中止操作。")
    print("!"*20)
    time.sleep(3)

    try:
        with Pool(processes=args.cores) as pool:
            list(tqdm(pool.imap_unordered(delete_file, files_to_process), total=num_to_delete, desc="[3/3] 正在并行删除文件", unit="file"))
    except Exception as e:
        print(f"删除过程中发生严重错误: {e}")
        sys.exit(1)

    print(f"\n操作完成。共删除了 {num_to_delete} 个文件。")

if __name__ == '__main__':
    main()
