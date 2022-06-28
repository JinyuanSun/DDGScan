#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="DDGScan",
    version="0.1.0",
    author="Jinyuan Sun",
    author_email="jysun@im.ac.cn",
    description="A powerful package for in silico mutation",
    url="https://github.com/JinyuanSun/DDGScan",
    packages=find_packages(),
    scripts=['DDGScan', 'grape-fast.py', 'multimer_scan.py'],
    install_requires=['pandas',
                      'numpy',
                      'joblib',
                      'sklearn',
                      'seaborn',
                      'matplotlib',
                      'venn',
                      'logomaker',
                      'modeller',
                      'openmm',
                      'mdtraj',
    ]
)
