#!/usr/bin/env python

import chardet
import sys
# 文件所在目录

files = sys.argv[1:]
for file in files:
    with open(file,'rb') as f:
        content=f.read()
    coding = chardet.detect(content).get("encoding") #获取文本编码信息
    print(f"{file}编码方式为{coding}")
    with open(file, "r", encoding="gbk") as f1, open("new-"+file, "w", encoding="utf-8") as f2:
        for l in f1:
            f2.write(l)