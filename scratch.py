#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 18:55:47 2020

@author: Pranay S. Yadav
"""
import itertools as it
from array import array
from collections import Counter
from functools import lru_cache

import numpy as np

import ETC
from ETC.nsrps_c import (
    N1,
    check_equality,
    get_entropy,
    get_mask_pairs,
    get_mask_windows,
    substitute_pairs,
    substitute_windows,
)

# from ETC.nsrws_core import get_mask_pairs as gmp
# from ETC.nsrws_core import substitute_pairs as sbp
# from ETC.nsrws_core import check_equality as ceq
# from ETC.nsrws_core import get_mask_windows as gmw
# from ETC.nsrws_core import substitute_windows as sbw


size = 20
np.random.seed(1)
y = np.random.randint(1, 3, size=size, dtype="u4")
idx = np.arange(1, len(y) - 1, 2, dtype="u4")
L = y.tolist()
z = array("I", y)
zC = array("I", [0] * (np.max(z) + 1))
q = array("I", tuple(range(100_000)))
order = 2

# def get_strided(y):
#     mask = masker(y)[:-1]
#     y_windowed = as_strided(y, shape=(len(y)-2+1,2), strides=y.strides*2, writeable=False)[mask,:]
#     uniques = np.unique(y_windowed, axis=0, return_counts=True, return_inverse=True)
#     idxes = np.where(uniques[1] == uniques[2].argmax())[0].astype('u4')
#     a = 1+y.max()
#     return substitute2(y, idxes, a)

# def get_stridedx(y):
#     mask = masker(y)[:-1]
#     y_windowed = as_strided(y, shape=(len(y)-2+1,2), strides=y.strides*2, writeable=False)[mask,:]
#     freq = Counter(tuple(x) for x in y_windowed.tolist()).most_common(1)[0]
#     pair = np.array(freq[0], dtype='u4')
#     a = 1+y.max()
#     return substitute3(y, pair, a)

# def get_stridedz(y):
#     mask = masker(y)[:-1]
#     y_windowed = as_strided(y, shape=(len(y)-2+1,2), strides=y.strides*2, writeable=False)[mask,:]
#     freq = Counter(tuple(x) for x in y_windowed.tolist()).most_common(1)[0]
#     pair = np.array(freq[0], dtype='u4')
#     a = 1+y.max()
#     return substitute4(y.copy(), pair, a)

# def X(y):
#     mask = xmasker(y)[:-1]
#     # y_windowed = as_strided(y, shape=(len(y)-2+1,2), strides=y.strides*2, writeable=False)[mask,:]
#     y_windowed = it.compress(zip(y[:-1], y[1:]), mask)
#     freq = Counter(y_windowed).most_common(1)[0]
#     pair = np.array(freq[0], dtype='u4')
#     a = 1+y.max()
#     xsubst(y, pair, a)
#     return y[y>0]

# def A(z):
#     mask = get_mask_pairs(z)[:-1]
#     z_windowed = it.compress(zip(z[:-1], z[1:]), mask)
#     freq = Counter(z_windowed).most_common(1)[0]
#     pair = array('I',freq[0])
#     a = 1+max(z)
#     substitute_pairs(z, pair, a)
#     return array('I', (x for x in z if x))


def nsrps(z):

    mask = get_mask_pairs(z)[:-1]
    z_windowed = it.compress(zip(z[:-1], z[1:]), mask)
    freq = Counter(z_windowed).most_common(1)[0]
    pair = array("I", freq[0])
    a = 1 + max(z)

    return array("I", substitute_pairs(z, pair, a))


def nsrps2(z):

    mask = gmp(z)[:-1]
    z_windowed = it.compress(zip(z[:-1], z[1:]), mask)
    freq = Counter(z_windowed).most_common(1)[0]
    pair = array("I", freq[0])
    a = 1 + max(z)

    return array("I", sbp(z, pair, a))


def xsrps(z, order):

    mask = get_mask_windows(z, order)[: -(order - 1)]
    z_windowed = it.compress(zip(*(it.islice(z, i, None) for i in range(order))), mask)
    freq = Counter(z_windowed).most_common(1)[0]
    window = array("I", freq[0])
    a = 1 + max(z)

    return array("I", substitute_windows(z, order, window, a))


def etc(x):
    z = x[:]
    etc = 0
    while not check_equality(z):

        z = nsrps(z)
        etc += 1
        # print(z)
    return etc


def etc2(x):
    z = x[:]
    etc = 0
    while not ceq(z):

        z = nsrps2(z)
        etc += 1
        # print(z)
    return etc


def xtc(x, order):
    z = x[:]
    etc = 0
    while not ceq(z) and len(z) >= order:

        z = xsrps(z, order)
        etc += 1
        # print(z)
    return etc


# %%

# np.random.seed(1)
# for n in [10_000, 50_000, 100_000]:
#     for m in range(50):
#         y = np.random.randint(1,3, size=n, dtype='u4')
#         z1 = array('I', y)
#         z2 = array('I', y)
#         z3 = tuple(y)
#         q1 = list(nsrps(z1))
#         q2 = list(xsrps(z2, 2))
#         q3 = ETC.NSRWS.run_once_NSRWS(z3, 2, 0)

#         print(n,m,q1==q2==q3, len(q1), len(q2), len(q3))
#         if not q1==q2==q3:
#             break
#     if not q1==q2==q3:
#             break


# %%


# results = []

# for n in [10, 100, 1000, 10_000, 50_000, 100_000, 200_000]:

#     np.random.seed(1)
#     y = np.random.randint(1,3, size=n, dtype='u4')
#     L = y.tolist()
#     z = array('I',y)

#     x = %timeit -o -r 10 -n 10 etc(z)
#     results.append({'method':'new', 'mean':x.average, 'sd':x.stdev, 'best':x.best, 'worst':x.worst, 'size':n})

#     x = %timeit -o -r 10 -n 10 ETC.compute(z, 2, False)
#     results.append({'method':'old', 'mean':x.average, 'sd':x.stdev, 'best':x.best, 'worst':x.worst, 'size':n})


# %%

case3 = (
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
)

case4 = (
    (1, 1, 1),
    (1, 1, 1),
    (1, 1, 2),
    (1, 2, 1),
    (2, 1, 1),
    (1, 1, 1),
    (1, 1, 1),
)


def test(win):

    idx = 2
    breaker = False
    skip = 0
    counter = Counter(dict(zip(set(win), it.repeat(0))))
    while idx <= len(win):

        print(idx - 2, idx - 1, end="\t")

        if not breaker:
            counter[win[idx - 2]] += 1
            if win[idx - 2] == win[idx - 1]:
                print(True, end="\t")
                if idx < len(win):
                    if win[idx - 2] == win[idx]:
                        skip = 2
                        print(True, skip)
                    else:
                        skip = 1
                        print(False, skip)
            else:
                skip = 0
                print(False, end="\t")
                if idx < len(win):
                    if win[idx - 2] == win[idx]:
                        breaker = True
                        print(True, breaker)
                    else:
                        breaker = False
                        print(False, breaker)

            # print(False, skip)

        idx += 1 + skip
    return counter


# @lru_cache(maxsize=None)
def compare(windows):
    """
    Returns a Boolean Mask where:
        First element is always True
        Remaining elements True if != first element, else False

    """
    out = [True if x != windows[0] else False for x in windows]
    out[0] = True

    return out


# @lru_cache(maxsize=None)
def slide(win, order=2):
    """
    Slides across a tuple of windows (tuples) and counts unique ones

    win : tuple of tuples of ints, each tuple element has length = order
    order : length of window to examine

    example input A for order=3:
        ((1,1,2), (1,2,1), (2,1,1))
        should return counter with 1 count for each window, & mask = [True, True, True]

    example input B for order=3:
        ((1,2,1), (2,1,2), (1,2,1))  # first and last are same
        should return 1 count each for 1st 2 windows, & mask = [True, True, False]

    example input C for order=3:
        ((1,1,1), (1,1,1), (1,1,1))  # all 3 are same
        should return counter with 1 count for 1st window, & mask = [True, False, False]
    """
    # initialize counter for unique tuples
    counter = Counter(dict(zip(set(win), it.repeat(0))))

    # initialize a boolean mask with True values
    mask = list(it.repeat(True, len(win)))

    # iterate over each window
    for n, window in enumerate(win):  # 1st loop
        print(n, window, mask[n], end="\t")
        # proceed only if mask is True for that window
        if mask[n]:

            # count that window
            counter[window] += 1

            # check if any successive windows are similar
            comp = compare(win[n : n + order])  # 2nd loop
            print(True, comp)
            # if any similar windows found, mask them
            if not all(comp):
                print("ALL")
                mask[n : n + order] = comp
        print()
    return counter, mask


def xlide2(win, order=2):

    # initialize counter for unique tuples
    counter = Counter(dict(zip(set(win), it.repeat(0))))

    # initialize a boolean mask with True values
    mask = list(it.repeat(True, len(win)))

    # iterate over each window
    for n, window in enumerate(win):  # 1st loop

        # proceed only if mask is True for that window
        if mask[n]:

            # count that window
            counter[window] += 1

            #
            mask[n + 1 : n + order] = [
                window != a and b
                for a, b in zip(win[n + 1 : n + order], mask[n + 1 : n + order])
            ]

    return counter, mask


def xompare(wind, mask):

    return [b and wind[0] != a for a, b in zip(wind[1:], mask[1:])]


@lru_cache(maxsize=None)
def clide(win, order=2):

    # initialize counter for unique tuples
    counter = Counter(dict(zip(set(win), it.repeat(0))))

    # initialize a boolean mask with True values
    mask = list(it.repeat(True, len(win)))

    # iterate over each window
    for n, window in enumerate(win):  # 1st loop

        # proceed only if mask is True for that window
        if mask[n]:

            # count that window
            counter[window] += 1

            #
            mask[n + 1 : n + order] = [
                window != a and b
                for a, b in zip(win[n + 1 : n + order], mask[n + 1 : n + order])
            ]

    return counter, mask


# %%
# order = 3x
# for n in range(6-order):    # 1st loop

#         # proceed only if mask is True for that window

#     for k in range(1,order):
#         for m in range(order):
#             print(n,k,m,'\t', n+m, n+k+m)
#         print()
#     print('-----')
