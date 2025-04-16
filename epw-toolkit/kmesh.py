#!/usr/bin/env python3
import sys

# 获取命令行参数
numargs = len(sys.argv) - 1

# 参数检查
if numargs < 3 or numargs > 4:
    print("usage: n1 n2 n3 [wan]")
    print("       n1  - divisions along 1st recip vector")
    print("       n2  - divisions along 2nd recip vector")
    print("       n3  - divisions along 3rd recip vector")
    print("       wan - omit the kpoint weight (optional)")
    sys.exit()

# 获取输入的 n1, n2, n3
n1, n2, n3 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

# 参数检查：确保 n1, n2, n3 都大于 0
if n1 <= 0:
    print("n1 must be > 0")
    sys.exit()
if n2 <= 0:
    print("n2 must be > 0")
    sys.exit()
if n3 <= 0:
    print("n3 must be > 0")
    sys.exit()

# 计算总的 k 点数量
totpts = n1 * n2 * n3

# 如果只有 3 个参数（不包括 wan）
if numargs == 3:
    print("K_POINTS crystal")
    print(totpts)
    for x in range(n1):
        for y in range(n2):
            for z in range(n3):
                # 格式化输出 k 点信息
                print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}{1/totpts:14.6e}")

# 如果有 4 个参数（包括 wan）
elif numargs == 4:
    for x in range(n1):
        for y in range(n2):
            for z in range(n3):
                # 格式化输出 k 点信息（没有权重）
                print(f"{x/n1:12.8f}{y/n2:12.8f}{z/n3:12.8f}")
