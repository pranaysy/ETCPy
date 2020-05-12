#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
cimport cython
import numpy as np
from libc.math cimport log2
cimport numpy as np

cpdef double get_entropy(unsigned int[::1] x):

    # cdef np.ndarray[np.npy_int64, ndim=1] counts = np.bincount(x)
    cdef long int[:] counts_view = np.bincount(x)
    # cdef long int[:] counts_view = counts
    cdef double counts_total = 0
    cdef Py_ssize_t counts_size = counts_view.shape[0]
    cdef Py_ssize_t m

    for m in range(counts_size):
        counts_total += counts_view[m]

    cdef double E = 0.0
    cdef double prob, logprob

    m = 0
    for n in range(counts_size):
        if counts_view[n]!=0:
            prob = counts_view[n] / counts_total
            logprob = log2(prob)
            E = E-prob*logprob

    return E
