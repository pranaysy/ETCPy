#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(
        ["./NSRWS/x1D/core.pyx", "./NSRWS/x2D/core.pyx", "./seq/estimates.pyx",],
        annotate=True,
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
    name="ETCPy",
    version="1.0",
    author_email="mail@pranaysy.com",
    description="Compute the Effort-To-Compress (ETC) of a symbolic sequence",
    packages=["ETC",],
    license="Apache License, Version 2.0",
    long_description=open("README.md").read()
)
