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
        [
            "ETC/NSRWS/x1D/compute_core.pyx",
            "ETC/NSRWS/x2D/compute_core.pyx",
            "ETC/common/compute_estimates.pyx",
        ],
        annotate=True,
        compiler_directives={"language_level": "3"},
    ),
    include_dirs=[numpy.get_include()],
)