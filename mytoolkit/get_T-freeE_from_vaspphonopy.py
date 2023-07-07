#!/usr/bin/env python
import os

# -o选项表示只输出匹配到的部分，而不是整行。这意味着grep命令只会输出符合正则表达式的部分，而不包括整行的其他内容。在我们的例子中，它将只输出匹配到的温度或自由能值，而不会包含"temperature:"或"free_energy:"等其他文本。
# -P选项表示使用Perl兼容的正则表达式（PCRE）进行匹配。Perl正则表达式支持更丰富的语法和功能，相比标准正则表达式，可以提供更强大的模式匹配能力。在我们的例子中，使用了一些Perl正则表达式的特性，如\K用于忽略前面的匹配内容，以及(\.\d+)?用于匹配可选的小数部分。
# \K 表示忽略前面的匹配内容
# ?是一个特殊字符，用于指示前面的元素是可选的。
# 匹配一个可选的负号（-?）
# (\.\d+)?用于匹配可选的小数部分。

T=os.popen("grep -oP 'temperature:\s+\K\d+' thermal_properties.yaml").read().split()
freeE=os.popen("grep -oP 'free_energy:\s+\K-?\d+(\.\d+)?' thermal_properties.yaml").read().split()
#[@]是用于引用整个数组array的所有元素
# #是获取数组长度的操作符。

with open("T-freeE.dat", "w") as tf:
    for t, f in zip(T, freeE):
        tf.write("{:>5} {:>12.7f}\n".format(int(t), float(f)*0.0103642696554971))

