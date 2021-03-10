#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(
        [
            "./ETC/NSRWS/x1D/core.pyx",
            "./ETC/NSRWS/x2D/core.pyx",
            "./ETC/seq/estimates.pyx",
            "./ETC/LZ76/core.pyx",
        ],
        annotate=False,
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
    name="ETCPy",
    version="1.3.5",
    author_email="mail@pranaysy.com",
    description="Compute the Effort-To-Compress (ETC) of a symbolic sequence",
    packages=find_packages(),
    license="Apache License, Version 2.0",
)
