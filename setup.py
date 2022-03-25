#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="GRAPE",
    version="1.0",
    author="Jinyuan Sun",
    author_email="jysun@im.ac.cn",
    description="A powerful package for GRAPE",
    url="https://github.com/JinyuanSun/DDGScan",
    packages=find_packages(),
    scripts=['grape-fast.py', 'gluster.py'],
    install_requires=['pandas',
                      'numpy',
                      'joblib',
                      'sklearn'
    ]
)
