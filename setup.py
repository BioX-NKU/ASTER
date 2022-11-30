#!/usr/bin/env python
#-*- coding:utf-8 -*-


from setuptools import setup, find_packages

setup(
    name="epiaster",
    version="0.0.5",
    keywords=("pip", "aster", "single-cell"),
    description="ASTER: accurate estimation of cell-type numbers in single-cell chromatin accessibility data",
    long_description="ASTER provides an accurate and efficient way to estimate the number of cell types in single-cell chromatin accessibility data. We provide documentation in the form of functional application programming interface documentation, tutorials and example workflows at https://aster.readthedocs.io/en/latest/index.html. All ASTER wheels distributed on PyPI are MIT licensed.",
    license="MIT Licence",
    url="https://github.com/BioX-NKU/ASTER",
    author="BioX-NKU",
    packages=find_packages(),
    python_requires='>=3.10.0',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        'numpy>=1.21.6,<1.22',
        'scipy>=1.8.0',
        'scikit-learn>=1.1.0',
        'kneed>=0.7.0',
        'scanpy>=1.9.1',
        'episcanpy==0.3.2',
        'igraph>=0.9.10',
        'louvain>=0.7.1',
        'leidenalg>=0.8.10',
        'tqdm'
    ]
)