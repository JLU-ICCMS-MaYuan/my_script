from setuptools import setup, find_packages

setup(
    name="my_script",
    version="2.0",
    author="madegan",
    author_email="myth620137018@163.com",
    description="一个计算qe vasp epw 和产生结构的小软件",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://gitee.com/mayuan_JLUPHY/my_script",
    project_urls={
        "Documentation": "https://your-docs-url",
        "Source": "https://gitee.com/mayuan_JLUPHY/my_script"
    },
    packages=find_packages(include=['qe', 'vasp', 'structuregenerator']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    install_requires=[
        "pymatgen==2022.9.21",
        "ase==3.25.0", # 我之前安装了3.24.1的ase的版本，无法读取relax.out文件，所以这里我强制要求安装3.22.1的版本
        "matplotlib",
        "f90wrap",
        "Cython",
        "scikit-learn",
        "pyfiglet==1.0.2",
        "spglib==2.5.0",
        "phonopy==2.28.0",
    ],
    scripts=[
        'qe/qe_main.py',
        'vasp/vasp_main.py',
        'structuregenerator/generator_main.py',
        'epw/epw_main.py',
    ],
    python_requires='>=3.7',
)