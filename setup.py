#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='guessmylt',
    version='0.2.2',

    description='An efficient way to guess the library type of RNA-reads',

    url='https://github.com/NBISweden/GUESSmyLT',
    author='Hampus Olin, Erik Berner-Wik, Caitlin Vigetun Haughey, Lisa Klasson, Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython==1.67', 'bcbio-gff==0.6.4', 'pysam>=0.13.0', 'snakemake==5.4.*'],
    include_package_data=True,

    entry_points={
        'console_scripts': ['GUESSmyLT = GUESSmyLT.GUESSmyLT:main',
        'GUESSmyLT-example = GUESSmyLT.test.example:main',
        ],
    },
)
