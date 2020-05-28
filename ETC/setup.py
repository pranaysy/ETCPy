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
        ["ETC/NSRWS/x1D/core.pyx", "ETC/NSRWS/x2D/core.pyx", "ETC/seq/estimates.pyx",],
        annotate=True,
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
)
