import os
from setuptools import setup

def parse_requirements():
    requires = []
    with open('requirements.txt', 'r') as fr:
        for line in fr:
            pkg = line.strip()
            if pkg.startswith('git+'):  # 如果是 git 链接，直接调用 pip 安装
                pip_install_git(pkg)
            else:
                requires.append(pkg)
    return requires

def pip_install_git(link):
    os.system(f'pip install --upgrade {link}')
    return

setup(
    name="my_script",
    version="2.0",
    author="madegan",
    author_email="myth620137018@163.com",
    description="一个计算qe vasp 和产生结构的小软件",
    url="https://gitee.com/mayuan_JLUPHY/my_script",
    packages=["qe", "vasp", "structuregenerator"],
    long_description="这个库有三个主要的功能, qe声子超导计算, vasp声子计算, 产生结构",
    classifiers=[
        'Development Status :: 1 - Alpha',
        'Intended Audience :: everybudy',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    install_requires=parse_requirements(),
    scripts=['qe/qe_main.py', 'vasp/vasp_main.py', 'structuregenerator/generator_main.py'],
    python_requires='>=3',
)
