from setuptools import setup, find_packages

setup(
    name="GRAPE",
    version="1.0",
    author="Jinyuan Sun",
    author_email="jysun@im.ac.cn",
    description="A powerful package for GRAPE",

    # 项目主页
    url="https://github.com/JinyuanSun/DDGScan", 

    # 你要安装的包，通过 setuptools.find_packages 找到当前目录下有哪些包
    packages=find_packages(),
    scripts=['grape-fast.py','gluster.py']
)
