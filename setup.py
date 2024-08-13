import os
import sys
from pathlib import Path
from setuptools import setup

def parse_requirements():
    requires = []
    with open('requirements.txt', 'r') as fr :
        for line in fr :
            pkg = line.strip()
            if pkg.startswith('git+'): # Python startswith() 方法用于检查字符串是否是以指定子字符串开头，如果是则返回 True，否则返回 False。如果参数 beg 和 end 指定值，则在指定范围内检查。
                pip_install_git(pkg)
            else:
                requires.append(pkg)
    return requires

def pip_install_git(link):
    os.system('pip install --upgrade {}'.format(link))
    return

#def Write_Tobin(qebin_path, vaspbin_path, my_scriptrc_ini):
#    with open(my_scriptrc_ini, "r") as f:
#        content = f.read()
#
#    with open(qebin_path, "w") as qe:
#        qe.writelines(content)
#
#    with open(vaspbin_path, "w") as vasp:
#        vasp.write(content)
#
#my_scriptrc_ini = Path.home().joinpath(".my_scriptrc.py")
#qebin_path      = Path("qe").joinpath("qebin.py")
#vaspbin_path    = Path("vasp").joinpath("vaspbin.py")
#Write_Tobin(qebin_path, vaspbin_path, my_scriptrc_ini)

setup(
    name="my_script",
    version="1.0",
    author="madegan",
    author_email="myth620137018@163.com",
    description="一个计算qe vasp 和产生结构的小软件",
    # 项目主页
    url="https://gitee.com/mayuan_JLUPHY/my_script", 
    # 你要安装的包, 通过 setuptools.find_packages 找到当前目录下有哪些包, 也可以自己指定有哪些包要被打包
    packages=["qe", "vasp", "structuregenerator"],
    # 指定要打包的模块文件
    # py_modules=["symbolspg.py"],
    # 详细描述
    long_description="这个库有三个主要的功能, qe声子超导计算, vasp声子计算, 产生结构",
    # classifiers 参数说明包的分类信息
    classifiers = [
        # 发展时期,常见的如下
        #   1 - Alpha
        #   2 - Beta
        #   3 - Production/Stable
        'Development Status :: 1 - Alpha',
        # 开发的目标用户
        'Intended Audience :: everybudy',
        # 属于什么类型
        'Topic :: Software Development :: Build Tools',
        # 许可证信息
        'License :: OSI Approved :: MIT License',
        # 目标 Python 版本
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    # 安装时需要安装的依赖包
    install_requires= parse_requirements(),
    # 指定可执行脚本,安装时脚本会被安装到系统 PATH 路径下
    scripts = ['qe/qe_main.py', 'vasp/vasp_main.py', 'structuregenerator/generator_main.py', "mytoolkit/tool_main.py"],
    python_requires='>=3',
)

