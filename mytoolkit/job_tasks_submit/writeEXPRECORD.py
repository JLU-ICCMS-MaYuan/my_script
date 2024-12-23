#!/usr/bin/env python3 
import argparse
import csv
from datetime import datetime

def append_to_record(content):
    # 获取当前的日期和时间
    current_time = datetime.now()
    current_date = current_time.strftime("%Y-%m-%d")
    current_time_str = current_time.strftime("%H:%M:%S")
    
    # 打开文件并追加内容（以 CSV 格式写入）
    with open("EXP.RECORD.csv", "a", newline='') as file:
        writer = csv.writer(file)
        # 如果文件为空，则写入表头
        if file.tell() == 0:  # 检查文件是否为空
            writer.writerow(["Date", "Time", "Content"])
        # 写入当前日期、时间和内容
        writer.writerow([current_date, current_time_str, content])

def main():
    # 使用 argparse 解析命令行输入
    parser = argparse.ArgumentParser(description="Append content to EXP.RECORD with timestamp.")
    parser.add_argument("content", type=str, help="The content to append to the EXP.RECORD file.")
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用函数将内容和时间追加到文件
    append_to_record(args.content)

if __name__ == "__main__":
    main()
