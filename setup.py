#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
DESCRIPTION = "Compute the Effort-To-Compress (ETC) of symbolic sequences"
DISTNAME = 'ETCPy'
PLATFORMS = ['posix']
MAINTAINER = "Pranay S. Yadav"
MAINTAINER_EMAIL = "mail@pranaysy.com"
URL = 'https://github.com/pranaysy/ETCPy/'
LICENSE = "Apache License, Version 2.0"
DOWNLOAD_URL = 'https://github.com/pranaysy/ETCPy/'
VERSION = '1.3.5'

from setuptools import setup, find_packages
from Cython.Build import cythonize
from pathlib import Path
import numpy

cython_modules = [*Path.cwd().rglob("*.pyx")]

setup(
    ext_modules=cythonize(
        cython_modules,
        annotate=False,
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
    name=DISTNAME,
    version=VERSION,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    platforms=PLATFORMS,
    url=URL,
    download_url=DOWNLOAD_URL,
    description=DESCRIPTION,
    long_description=DESCRIPTION,
    packages=find_packages(),
    license=LICENSE,
    install_requires = ["numpy", "cython"],
    classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 3',
              'License :: OSI Approved :: Apache Software License',
              'Topic :: Scientific/Engineering :: Mathematics',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
)
