#!/usr/bin/env python3
import re

# 文件路径
filename = 'split_ph.in'

# 读取文件内容
with open(filename, 'r') as file:
    content = file.readlines()

# 查找并修改trans的值
for i, line in enumerate(content):
    match = re.search(r'\btrans\s*=\s*(\.(true|false)\.)', line, re.IGNORECASE) #re.IGNORECASE：忽略大小写，这意味着模式中的true或false可以是大写或小写形式。
    # \b单词边界（word boundary），确保匹配的是一个完整的单词，而不是单词的一部分。
    # (\.(true|false)\.)：匹配.true.或.false.，括号()表示这是一个捕获组（capture group），可以提取匹配的内容。
    # \.：匹配一个字面的点号.
    # (true|false)：匹配字符串true或false，竖线|表示“或”的意思
    # \.：匹配另一个字面的点号.
    if match:
        current_value = match.group(1).lower()
        new_value = '.false.' if current_value == '.true.' else '.true.'
        content[i] = re.sub(r'\btrans\s*=\s*\.(true|false)\.', f'trans = {new_value}', line, flags=re.IGNORECASE)
        break

# 写回修改后的文件内容
with open(filename, 'w') as file:
    file.writelines(content)

print(f"{content[i]}", end='')
