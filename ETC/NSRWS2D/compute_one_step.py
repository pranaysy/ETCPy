#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from array import array
from collections import Counter
from itertools import compress, islice
from time import perf_counter

from ETC.NSRWS2D import compute_core as cc


def _one_step_pairs(seq_x, seq_y, verbose=True):
    before = perf_counter()
    mask = cc.get_mask_pairs(seq_x, seq_y)
    seq_pairs_filt = compress(
        zip(zip(seq_x[:-1], seq_y[:-1]), zip(seq_x[1:], seq_y[1:])), mask
    )
    freq_pair, count = Counter(seq_pairs_filt).most_common(1)[0]
    sub_value_x = 1 + max(seq_x)
    sub_value_y = 1 + max(seq_y)
    pair_x = array("I", [freq_pair[0][0], freq_pair[1][0]])
    pair_y = array("I", [freq_pair[0][1], freq_pair[1][1]])
    if count == 1:
        out_x = array("I", seq_x[1:])
        out_x[0] = sub_value_x
        out_y = array("I", seq_y[1:])
        out_y[0] = sub_value_y
        signal = True
    else:
        out_x, out_y = cc.substitute_pairs(
            seq_x, seq_y, pair_x, pair_y, sub_value_x, sub_value_y
        )
        out_x = array("I", out_x)
        out_y = array("I", out_y)
        signal = False
    after = perf_counter()
    if verbose:
        return out_x, out_y, pair_x, pair_y, count, after - before, signal
    return out_x, out_y, signal


# def _one_step_windows(seq, order, verbose=True):

#     before = perf_counter()
#     mask = cc.get_mask_windows(seq, order)[: -(order - 1)]
#     z_windowed = compress(zip(*(islice(seq, i, None) for i in range(order))), mask)
#     z_windowed = tuple(z_windowed)
#     freq_window, count = Counter(z_windowed).most_common(1)[0]
#     sub_value = 1 + max(seq)
#     window = array("I", freq_window)
#     if count == 1:
#         out = array("I", seq[order - 1 :])
#         out[0] = sub_value
#         signal = True
#     else:
#         out = array("I", cc.substitute_windows(seq, order, window, sub_value))
#         signal = False
#     out = array("I", cc.substitute_windows(seq, order, window, sub_value))
#     after = perf_counter()
#     if verbose:
#         return out, freq_window, count, after - before, signal
#     return out, signal


def _one_step_windows(seq, order, verbose=True):
    pass


def _one_step(seq_x, seq_y, order, verbose=True):

    if order == 2:
        return _one_step_pairs(seq_x[:], seq_y[:], verbose)
    if order > 2:
        return _one_step_windows(seq_x[:], seq_y[:], order, verbose)


def one_step(seq_x, seq_y, order, verbose=True, check=True):
    if not isinstance(seq_x, array):
        seq_x = array("I", seq_x)
    if not isinstance(seq_y, array):
        seq_y = array("I", seq_y)

    if check and cc.check_equality(seq_x):
        print("> All elements in sequence x are equal!")
        return None
    if check and cc.check_equality(seq_y):
        print("> All elements in sequence y are equal!")
        return None

    if len(seq_x) < order or len(seq_y) < order:
        print("> Sequence input shorter than order!\n> Can't perform substitution ...")
        return None

    return _one_step(seq_x, seq_y, order, verbose)
