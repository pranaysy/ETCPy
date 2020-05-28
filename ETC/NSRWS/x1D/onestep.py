#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from array import array
from collections import Counter
from itertools import compress, islice
from time import perf_counter

from ETC.NSRWS.x1D import core as cc

def _apply_mask(seq, mask, order):

    filtered = compress(zip(*(islice(seq, i, None) for i in range(order))), mask)
    freq_window, count = Counter(filtered).most_common(1)[0]

    return freq_window, count

def _onestep_pairs(seq, verbose=True):
    before = perf_counter()
    mask = cc.get_mask_pairs(seq)
    freq_pair, count = _apply_mask(seq, mask, 2)
    sub_value = 1 + max(seq)
    pair = array("I", freq_pair)
    if count == 1:
        out = array("I", seq[1:])
        out[0] = sub_value
        signal = True
    else:
        out = array("I", cc.substitute_pairs(seq, pair, sub_value))
        signal = False
    after = perf_counter()
    if verbose:
        return out, freq_pair, count, after - before, signal
    return out, signal


def _onestep_windows(seq, order, verbose=True):

    before = perf_counter()
    mask = cc.get_mask_windows(seq, order)
    freq_window, count = _apply_mask(seq, mask, order)
    sub_value = 1 + max(seq)
    window = array("I", freq_window)
    if count == 1:
        out = array("I", seq[order - 1 :])
        out[0] = sub_value
        signal = True
    else:
        out = array("I", cc.substitute_windows(seq, order, window, sub_value))
        signal = False
    out = array("I", cc.substitute_windows(seq, order, window, sub_value))
    after = perf_counter()
    if verbose:
        return out, freq_window, count, after - before, signal
    return out, signal


def _onestep(seq, order, verbose=True):
    if order == 2:
        return _onestep_pairs(seq[:], verbose)
    if order > 2:
        return _onestep_windows(seq[:], order, verbose)


def onestep(seq, order, verbose=True, check=True):
    if not isinstance(seq, array):
        seq = array("I", seq)
    if check and cc.check_equality(seq):
        print("> All elements in sequence are equal!")
        return None
    if len(seq) < order:
        print("> Sequence input shorter than order!\n> Can't perform substitution ...")
        return None

    return _onestep(seq, order, verbose)
