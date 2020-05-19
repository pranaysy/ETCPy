#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

# import ETC
from ETC import compute_1D, compute_2D
from ETC.utils import partition
from functools import partial

get1D = partial(compute_1D, order=2, verbose=False, truncate=True)
get2D = partial(compute_2D, order=2, verbose=False, truncate=True)

size = 1000
# seq_x = ETC.generate(size=size)
# seq_y = ETC.generate(size=size)

# Past window size
LEN_past = 150

# Current window size
ADD_meas = 15

# Jump step size
STEP_size = 20

# Number of bins
Num_bins = 2

# Length of sequence
LEN = size


def compute(seq_x, seq_y, LEN_past, ADD_meas, STEP_size, partitions=False):

    if partitions:
        seq_x = partition(seq_x, partitions)
        seq_y = partition(seq_y, partitions)

    l_1D = []
    l_2D = []
    LEN = len(seq_x)
    LEN_to_check = LEN_past + ADD_meas

    for k in range(0, LEN - LEN_to_check, STEP_size):

        ETC1D_ini = compute_1D(seq_x[k : k + LEN_past], 2, 0, 1)["ETC1D"] / (
            LEN_past - 1
        )
        # print(ETC1D_ini, end=' ')

        ETC2D_ini = compute_2D(
            seq_x[k : k + LEN_past], seq_y[k : k + LEN_past], 2, 0, 1
        )["ETC2D"] / (LEN_past - 1)
        # print(ETC2D_ini, end=' ')

        ETC1D_fin = compute_1D(seq_x[k : k + LEN_to_check], 2, 0, 1)["ETC1D"] / (
            LEN_to_check - 1
        )
        # print(ETC1D_fin, end=' ')

        ETC2D_fin = compute_2D(
            seq_x[k : k + LEN_to_check],
            seq_y[k : k + LEN_past] + seq_x[k + LEN_past : k + LEN_to_check],
            2,
            0,
            1,
        )["ETC2D"] / (LEN_to_check - 1)
        # print(ETC2D_fin, end='\n\n')

        ETC1D_delta = ETC1D_fin - ETC1D_ini
        ETC2D_delta = ETC2D_fin - ETC2D_ini

        l_1D.append(ETC1D_delta)
        l_2D.append(ETC2D_delta)

    CCC = (sum(l_1D) - sum(l_2D)) / (len(l_1D) - 1)
    print(CCC)
    return CCC
