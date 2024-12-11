#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import io  # 导入io模块

# 使用argparse解析命令行参数
def parse_args():
    parser = argparse.ArgumentParser(description="读取log.lammps文件并绘制Time与TotEng的图")
    parser.add_argument('-l', '--log-file', type=str, help='log.lammps文件路径')
    parser.add_argument('-b', '--begin-line', type=int, default=0, help='跳过前面的行数，默认为0')
    parser.add_argument('-e', '--end-line', type=int, default=0, help='跳过后面的行数，默认为0')
    parser.add_argument('-o', '--output-plot', type=str, default='EngTot.png', help='保存图像为png文件的路径，默认为"EngTot.png"')
    parser.add_argument('-f', '--output-file', type=str, default='traj.csv', help='保存DataFrame数据为csv文件的路径，默认为"traj.csv"')
    parser.add_argument('-r', '--read-csv', type=str, help='读取CSV文件并绘制图，文件路径')
    return parser.parse_args()

# 读取CSV文件
def read_csv_file(csv_file):
    df = pd.read_csv(csv_file)
    print(df)
    return df

# 读取log.lammps文件
def read_log_file(log_file, begin_line, end_line):
    # 打开文件并读取所有行
    with open(log_file, 'r') as file:
        lines = file.readlines()
    
    # 跳过后面的skip_behindlines行
    lines = lines[begin_line-1:end_line] if (begin_line > 0 and end_line > 0) else lines

    # 通过pandas读取跳过行后的数据
    df = pd.read_csv(io.StringIO("".join(lines)), sep=r'\s+', header=None, 
                     names=['Step', 'Time', 'Temp', 'Press', 'PotEng', 'KinEng', 'TotEng', 'Lx', 'Ly', 'Lz'])
    #print(df.head(), "\n...\n...\n", df.tail().to_string(header=False))    
    print(df)    
    return df

# 绘制Time与TotEng的图
def plot_time_vs_toteng(df, output_plot=None):
    plt.plot(df['Time'], df['PotEng'], label='PotEng', color='r')
    #plt.plot(df['Time'], df['KinEng'], label='KinEng', color='g')
    plt.plot(df['Time'], df['TotEng'], label='TotEng', color='b')
    
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.title('Time vs Energies (PotEng, KinEng, TotEng)')
    plt.legend()  # 显示图例
    plt.grid(True)
    
    # 如果指定了输出路径，就保存为png文件
    if output_plot:
        plt.savefig(output_plot)
        print(f"图像已保存为 {output_plot}")
    else:
        plt.show()

# 保存DataFrame数据为csv文件
def save_data_to_csv(df, output_file=None):
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"数据已保存为 {output_file}")

def main():
    # 解析命令行参数
    args = parse_args()
    
    # 如果指定了读取CSV文件
    if args.read_csv:
        # 读取CSV文件并绘制图
        df = read_csv_file(args.read_csv)
        plot_time_vs_toteng(df, args.output_plot)
    else:
        # 否则，读取log文件
        df = read_log_file(args.log_file, args.begin_line, args.end_line)
        plot_time_vs_toteng(df, args.output_plot) 
        # 保存数据
        save_data_to_csv(df, args.output_file)

if __name__ == "__main__":
    main()
