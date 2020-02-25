#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: pranay
"""

import numpy as np
from numba import njit

@njit
def _discretize(x, n_bins):

    delta_inv = n_bins / (np.ptp(x) + 1e-6)

    return np.floor((x - x.min()) * delta_inv) + 1

def _partition(x, n_bins):
    
    return _discretize(x, n_bins).astype(np.int)

def partition(x, n_bins):
    
    if not isinstance(n_bins, int) or n_bins < 2:
        print('>> Number of bins is invalid ...')
        return None
    
    if check_inputs(x) or n_bins > len(x):
        return _partition(x, n_bins)
    else:
        return None
    
      
def generate(size=10, partitions=2):
    
    if not isinstance(partitions, int) or partitions < 2 or partitions > size:
        print('>> Number of bins is invalid ...')
        return None
    
    random = np.random.random(size)
    
    return _partition(random, partitions)

def warmup():
    
    from .NSRPS import one_step
    
    if one_step(generate()).any():
        print('>> Numba optimizations ready ...')
        
        if np.array_equal(np.array([3,3]),one_step(np.array([1,2,1,2]))):
            print('>> NSRPS ready ...')
            return True
        else:
            print('>> NSRPS failed ...')
            return False
    else:
        print('>> Numba optimizations failed ...')
        return False
    
def check_equality(x):
    
    if np.all(x == x[0]):
        return True
    else:
        return False
    
def check_inputs(x):
    
    if isinstance(x, (np.ndarray, np.number)):
        if len(x) > 1 and len(x.shape) == 1:
            print('>> Input has valid type and dimensions ...')
            return True
        else:
            print('>> Input has invalid dimensions ...')
            return False
    
    else:
        print('>> Input is not a numeric NumPy array ...')
        return False

def check_sanity(x):
    if check_inputs(x) and not check_equality(x):
        return True
    else:
        return False