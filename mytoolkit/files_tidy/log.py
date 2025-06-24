#!/usr/bin/env python3
import os
import sys

def main():
    # 获取参数
    if len(sys.argv) > 1:
        # 如果有参数，使用参数作为补充说明
        additon_info = sys.argv[1]
    else:
        additon_info = ''
        
    # 确定目标文件路径
    home_dir = os.path.expanduser('~')
    workdir_file = os.path.join(home_dir, 'mylog')
    
    # 追加路径到文件
    path_to_append = os.getcwd()
    try:
        with open(workdir_file, 'a') as f:
            f.write(path_to_append + '  ' + additon_info +'\n' )
        print(f"Successfully appended: {path_to_append + '  ' + additon_info}")
    except IOError as e:
        print(f"Error: Could not write to {workdir_file}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
