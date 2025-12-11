from setuptools import setup, find_packages

# 使用 pyproject.toml 作为主配置，setup.py 仅保留兼容入口，避免与 PEP 621 配置重复。
setup(
    name="my_script",
    version="2.0.0",
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
    packages=find_packages(include=['qe', 'vasp', 'epw', 'structuregenerator']),
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
    python_requires='>=3.7',
)
