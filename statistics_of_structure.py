#!/work/home/may/miniconda3/bin/python3
'''
Use examples:
    statistics_of_structure.py -d 195/ -o -hc 75 100 -p 75 100
    -d 195/                 结构所在的母目录,
                            这个目录里面必须包含一个名为struc的子目录
                            这个目录里面必须包含spap运行后的输出文件
    -o -hc 75 100           输出 氢含量在某一个范围内的所有结构 到一个名为7 5.0-100.0 的目录
    -p 75 100               输出 氢含量在某一个范围内的所有结构 的对称性统计图
'''

import os
import re
import shutil
import pandas as pd


from collections import defaultdict
from argparse import ArgumentParser


def read_spap_file(src):

    contents_ao = open(os.path.join(src, "Analysis_Output.dat" ), "r").readlines()
    contents_ss = open(os.path.join(src, "structure_source.dat"), "r").readlines()

    white_space  = r"\s+"
    space_group  = r"\w*\-*\w*[(]*\d+[)]"                                 ; spg_patter   = re.compile(space_group)
    number_index = r"\s*\d{1,5}"                                          ; ni_patter    = re.compile(number_index)
    stoichometry = r"[a-zA-Z]{1,2}\d+[a-zA-Z]{1,2}\d+"                    ; stoi_patter  = re.compile(stoichometry)
    file_name    = r"\d*\-*[a-zA-Z]{1,2}\d+[a-zA-Z]{1,2}\d+\-{3}\w*\.vasp"; fn_patter    = re.compile(file_name)


    struct_info_ss = defaultdict(dict)
    struct_info_ao = defaultdict(dict)
    for cont in contents_ss:
        ni           = re.search(ni_patter, cont).group().strip(" ")
        stoi_res     = re.search(stoi_patter, cont).group()
        file_name    = re.search(fn_patter, cont).group()
        X_num, H_num = re.findall(r"\d+", stoi_res)

        H_content    = round(int(H_num) / (int(X_num) + int(H_num)), ndigits=3)
        struct_info_ss[ni]["stoichiometry"] = stoi_res
        struct_info_ss[ni]["filename"]      = file_name
        struct_info_ss[ni]["H_content"]     = round(H_content * 100, ndigits=3)

    for cont in contents_ao:
        spg_res = re.search(spg_patter, cont)
        ni_res  = re.search(ni_patter, cont)
        if spg_res is not None and ni_res is not None:
            spg_name = spg_res.group().split("(")[0]
            spg_num  = spg_res.group().split("(")[1].strip(")").strip(" ")
            ni = re.search(ni_patter, cont).group().strip(" ")

            struct_info_ao[ni]["spg_name"]      = spg_name
            struct_info_ao[ni]["spg_num"]       = int(spg_num)
            struct_info_ao[ni]["stoichiometry"] = struct_info_ss[ni]["stoichiometry"]
            struct_info_ao[ni]["H_content"]     = struct_info_ss[ni]["H_content"]
            struct_info_ao[ni]["filename"]      = struct_info_ss[ni]["filename"]

    return struct_info_ao

def get_pandas(struct_info_ao):
    # spg_set = {value["spg_num"] for value in struct_info_ao.values()}

    df = pd.DataFrame(struct_info_ao).T
    # print(df)
    df.sort_values(
        by=["spg_num", "H_content"],
        ascending=[False, False],
        inplace=True
    )

    return df

def W_eccel_O_struct(df, src, h_low, h_high):
    if os.path.exists(os.path.join(src, "all_struct.xlsx")):
        os.remove(os.path.join(src, "all_struct.xlsx"))
    excelname = os.path.join(src, "all_struct_info.xlsx")
    df.to_excel(excelname, sheet_name=destination_dir.strip("/"))  # 这里是一个很容易疏忽的地方可能导致bug。
                                                                   # 如果 -d 201 sheet_name=destination_dir 就够用了
                                                                   # 如果 -d 201/ sheet_name=destination_dir.strip("/") 才够用

    if not h_high >= h_low:
        h_high, h_low = h_low, h_high
    output_path = os.path.join(src, str(h_low)+"-"+str(h_high))
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    for index, row in df.iterrows():
        if h_low <= row["H_content"] < h_high :
            src_file = os.path.join(sub_struc_dir, row["filename"])
            des_file = os.path.join(output_path, "get"+str(row["spg_num"])+"-"+"from"+row["filename"])
            try:
                shutil.copy(src_file, des_file)
            except:
                print(src_file, "doesn't exist! ")

def plot_p(df, src, h_low, h_high):
    import matplotlib.pyplot as plt

    # 按照H含量筛选, H_content在75%以上时，结构都有哪些
    if not h_low <= h_high:
        h_low, h_high = h_high, h_low
    df_hcont   = df[(df["H_content"] >= h_low) & (df["H_content"] < h_high)]
    # 再在此基础上统计，这些结构的对称性分布有什么特征
    df_spg_num = df_hcont["spg_num"].value_counts()
    # 设置画布大小
    plt.figure(dpi=150, figsize=(8,7)); 
    df_spg_num.plot(kind='bar', width = 0.2, title="Space group symmetry statistics for structures\n with hydrogen content between {}% and {}%".format(h_low, h_high))
    plt.xlabel("Space group number");   plt.xticks(rotation=-60)
    plt.ylabel("Quantity")          ;   plt.grid(axis="y", c='grey', linestyle='--', which="major")
    plt.savefig(os.path.join(src, "Hydrogen-content-distribution-diagram"))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--destination_dir",
        default=None,
        type=str,
        help="请输入结构所在的母目录, 改母目录中必须包含struc目录"
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store_true",
        default=False,
        help="是否输出所有不相似的结构"
    )
    parser.add_argument(
        '-hc',
        '--hydrogen-content',
        action='store',
        type=float,
        dest='hydrogen_content',
        # required=True,
        nargs='+',
        default=[None, None]
    )
    parser.add_argument(
        '-p',
        '--plot-picture',
        action='store',
        type=float,
        dest='plot_picture',
        # required=True,
        nargs='+',
        default=[75, 100],
        help="是否绘制统计分布图, 如果绘制, 输入范围"
    )
    args = parser.parse_args()
    destination_dir = args.destination_dir
    output          = args.output
    h_low, h_high   = args.hydrogen_content
    ph_low, ph_high = args.plot_picture

    src             = os.path.abspath(destination_dir)

    sub_struc_dir   = os.path.join(src, "struc")

    struct_info_ao  = read_spap_file(src)
    df              = get_pandas(struct_info_ao)
    if output and h_low and h_high:
        W_eccel_O_struct(df, src, h_low, h_high)
    
    if ph_low and ph_high:
        plot_p(df, src, ph_low, ph_high)

