#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: pranay
"""

import numpy as np
from numba import njit
from .utils import check_equality, check_sanity

def _get_all_pairs(x):
    
    y = np.vstack((x[:-1], x[1:]))
    M = np.int((y.max() + 1))
    
    uniqs, invind, counts = np.unique(M*y[0] + y[1], return_counts=True, return_inverse=True)
       
    matr = np.zeros(M*M)
    matr[uniqs] = counts
    uniqc = matr.reshape(M,M)

    return uniqs, invind, uniqc

@njit
def _get_freq_pair_idx(x):
    
    uniqs, invind, uniqc = x
    np.fill_diagonal(uniqc, 0)
    
    return np.where(invind == np.where(uniqs == uniqc.argmax())[0])[0]

@njit
def _substitute(x, idx):
    
    y = x.copy()
    y[idx] = y.max() + 1
    filtidx = np.full(y.shape, True)
    filtidx[idx+1] = False
    
    return y[filtidx]

def one_step(x):
    
    return _substitute(x, _get_freq_pair_idx((_get_all_pairs(x))))

def NSRPS(x):
    
    if not check_sanity(x):
        return False
    
    if check_equality(x):
        print('> Constant symbolic sequence, nothing to substitute ...')
        return None
    else:
        return one_step(x)