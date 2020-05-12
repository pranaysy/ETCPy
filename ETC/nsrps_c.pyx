#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
import array
import itertools as it
from collections import Counter

import numpy as np

from cpythoncimportarray import bool

cimport numpy as np
cimport cython

# DTYPE = np.uintc
#  ctypedef np.uint32_t DTYPE_t

# @cython.boundscheck(False)
# @cython.wraparound(False)
# cpdef masker(np.ndarray[np.npy_uint32, ndim=1] x):

#     cdef np.ndarray[np.npy_bool, ndim=1] mask = np.ones(x.shape[0], dtype=np.bool)
#     cdef long int n = 0
#     cdef long int upperlim = x.shape[0] - 1

#     while n < upperlim:

#         if x[n] == x[n+1]:
#             if x[n+1] == x[n+2]:
#                 mask[n+1] = False
#                 n += 1
#         n += 1

#     return mask

# cpdef count(unsigned int[::1]x):
#     cdef Py_ssize_t x_max = len(x) - 1
#     cdef unsigned int a
#     cdef unsigned int b
#     cdef unsigned int xm = max(x) +1
#     cdef list counts = [[0 for _ in range(xm)] for _ in range(xm)]

#     for n in range(x_max):
#         a = x[n]
#         b = x[n+1]
#         counts[a][b] += 1

#     return counts


#



# @cython.boundscheck(False)
# @cython.wraparound(False)
# def substitute2(np.ndarray[np.npy_uint32, ndim=1] x, np.ndarray[np.npy_uint32, ndim=1] idx, long int a):

#     x[idx] = a
#     x[idx+1] = 0

#     return x

# @cython.boundscheck(False)
# @cython.wraparound(False)
# cpdef substitute4(np.ndarray[np.npy_uint32, ndim=1] x, np.ndarray[np.npy_uint32, ndim=1] pair, long int a):

#     cdef size_t n
#     cdef size_t x_size = x.shape[0]
#     for n in range(x_size):
#         if x[n] == pair[0] and x[n+1] == pair[1]:
#             x[n] = a
#             x[n+1] = 0

#     return x

# @cython.boundscheck(False)
# @cython.wraparound(False)
# cpdef substitute3(np.ndarray[np.npy_uint32, ndim=1] x, np.ndarray[np.npy_uint32, ndim=1] pair, long int a):
#     cdef np.ndarray[np.npy_uint32, ndim=1] out = x.copy()
#     cdef int n
#     cdef Py_ssize_t x_size = out.shape[0]
#     for n in range(x_size):
#         if out[n] == pair[0] and out[n+1] == pair[1]:
#             out[n] = a
#             out[n+1] = 0

#     return out


# @cython.boundscheck(False)
# @cython.nonecheck(False)
# @cython.wraparound(False)
# cpdef xmasker(unsigned int[::1] x):
#     cdef Py_ssize_t x_max = x.shape[0]
#     cdef np.ndarray[np.npy_bool, ndim=1] mask = np.ones(x_max, dtype=np.bool)
#     cdef long int n = 0
#     cdef long int upperlim = x.shape[0] - 1

#     while n < upperlim:

#         if x[n] == x[n+1]:
#             if x[n+1] == x[n+2]:
#                 mask[n+1] = False
#                 n += 1
#         n += 1

#     return mask

# @cython.boundscheck(False)
# @cython.nonecheck(False)
# @cython.wraparound(False)
# cpdef void xsubst(unsigned int[::1] x, unsigned int[::1] pair, long int a):
#     cdef Py_ssize_t n
#     cdef Py_ssize_t x_size = x.shape[0]

#     for n in range(x_size):
#         if x[n] == pair[0] and x[n+1] == pair[1]:
#             x[n] = a
#             x[n+1] = 0


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef array.array get_mask_pairs(unsigned int[::1] x):

    cdef Py_ssize_t x_size = len(x)

    cdef array.array int_template = array.array('I', [])
    cdef array.array mask = array.clone(int_template, x_size, zero=True)
    cdef unsigned int[:] mask_view = mask

    cdef Py_ssize_t n = 0
    cdef Py_ssize_t upperlim = x_size - 1

    for n in range(x_size):
        mask_view[n] += 1

    n = 0
    while n < upperlim:

        if x[n] == x[n+1] and x[n+1] == x[n+2]:

            mask_view[n+1] = 0
            n += 1

        n += 1

    return mask

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef list substitute_pairs(unsigned int[::1] x, unsigned int[::1] pair, long int a):

    cdef Py_ssize_t n
    cdef Py_ssize_t x_size = len(x)
    cdef list out = []

    for n in range(x_size-1):

        if x[n] == pair[0] and x[n+1] == pair[1]:

            x[n] = a
            x[n+1] = 0

    for n in range(x_size):

        if x[n]:

            out.append(x[n])

    return out

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
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

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef bool check_equality(unsigned int[::1] x):

    cdef Py_ssize_t n
    cdef Py_ssize_t x_size = len(x)
    cdef bool yes = True
    cdef bool no = False

    for n in range(x_size):
        if x[0] != x[n]:
            return no

    return yes



@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef array.array get_mask_windows(unsigned int[::1] x, unsigned int order):

    cdef Py_ssize_t x_size = len(x)

    cdef array.array int_template = array.array('I', [])
    cdef array.array mask = array.clone(int_template, x_size, zero=True)
    cdef unsigned int[:] mask_view = mask

    cdef Py_ssize_t n = 0
    cdef Py_ssize_t m = 0
    cdef Py_ssize_t k = 0
    cdef unsigned int track = 0
    cdef Py_ssize_t upperlim = x_size - 1

    for n in range(x_size):
        mask_view[n] += 1

    # iterate over each window
    for n in range(x_size-order):    # 1st loop

        # proceed only if mask is True for that window
        if mask_view[n]:
            for k in range(1,order):
                for m in range(order):
                    if x[n+m] == x[n+m+k]:
                        track += 1
                        # print(n+m,x[n+m], n+m+k,x[n+m+k], track)
                mask_view[n+k] = track!=order and mask_view[n+k]
                track = 0
    return mask

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cpdef list substitute_windows(unsigned int[::1] x, unsigned int order, unsigned int[::1] window, long int a):

    cdef Py_ssize_t n, m
    cdef unsigned int track = 0
    cdef Py_ssize_t x_size = len(x)
    cdef list out = []

    for n in range(x_size-order+1):

        for m in range(order):

            if x[n+m] == window[m]:
                track += 1

        if track == order:
            x[n] = a
            for m in range(1, order):
                x[n+m] = 0

        track = 0
    n = 0
    for n in range(x_size):

        if x[n]:

            out.append(x[n])

    return out

def N1(z):

    mask = get_mask_pairs(z)[:-1]
    z_windowed = it.compress(zip(z[:-1], z[1:]), mask)
    freq = Counter(z_windowed).most_common(1)[0]
    pair = array.array('I',freq[0])
    a = 1+max(z)
    out = array.array('I', substitute_pairs(z, pair, a))

    return out
