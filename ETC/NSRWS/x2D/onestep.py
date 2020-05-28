#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter
from itertools import compress, islice
from time import perf_counter

from ETC.NSRWS.x2D import core as cc
from ETC.seq.recode import cast
from ETC.seq.check import arraytype

def _mask_and_count(seq_x, seq_y, mask, order):
    """
    Apply binary mask to a pair of sequences & count most frequently jointly occurring windows

    This function does 3 things in the following sequence:
        1. Create sliding windows of a given size (order) - using zip and islice
        2. Apply a supplied mask to the sliding windows - using compress
        3. Count most frequently occurring window - using Counter

    In the NSRWS algorithm, this is the most time consuming step. Essentially expands
    two 1D sequences to a 2D sequence - where the sequences follows along row-wise &
    the columnar expansion encodes a sliding window for each row jointly from both
    sequences:
        1D sequences:
            (1,2,3,4,5,6,7)
            (3,4,5,6,7,8,9)

        2D expansion for window order=3:
            (((1,3),(2,4),(3,5)),
             ((2,4),(3,5),(4,6)),
             ((3,5),(4,6),(5,7)),
             ((4,6),(5,7),(6,8)),
             ((5,7),(6,8),(7,9)))

        The mask is applied row-wise & must be of the same length as the number of rows
        in this 2D expansion. This is given by:
            len(mask) = len(seq) - (order - 1)

        Example application of the mask (1,0,0,1,1):
            1 -> (((1,3),(2,4),(3,5)),
            0 ->  ((2,4),(3,5),(4,6)),          (((1,3),(2,4),(3,5)),
            0 ->  ((3,5),(4,6),(5,7)),    --->   ((4,6),(5,7),(6,8)),
            1 ->  ((4,6),(5,7),(6,8)),           ((5,7),(6,8),(7,9)))
            1 ->  ((5,7),(6,8),(7,9)))

        Unique windows (rows of 2D expansion) are counted and most frequently occurring
        row is returned with counts.

    Parameters
    ----------
    seq_x : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    seq_y : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    mask : array.array
        Collection of Booleans, where 0s indicate locations on "seq" to mask out.
        0s correspond to overlapping windows.
    order : int
        Size of window for NSRWS, 2 or greater.

    Returns
    -------
    pair_x : array.array
        Most frequently occurring non-overlapping "window" of size "order" in seq_x
    pair_y : array.array
        Most frequently occurring non-overlapping "window" of size "order" in seq_y
    count : int
        Number of times the most frequently occurring window occurs.

    """
    # Create overlapped sliding windows (each window a tuple of size order) & apply mask
    filtered = compress(zip(*(islice(zip(seq_x, seq_y), i, None) for i in range(order))), mask)

    # Count sliding windows (tuples are hashable!) & get the one most common with counts
    freq_pair, count = Counter(filtered).most_common(1)[0]

    # Assign array type and return
    pair_x = cast([freq_pair[0][0], freq_pair[1][0]])
    pair_y = cast([freq_pair[0][1], freq_pair[1][1]])

    return pair_x, pair_y, count

def _onestep_pairs(seq_x, seq_y, verbose=True):
    before = perf_counter()
    mask = cc.get_mask_pairs(seq_x, seq_y)

    pair_x, pair_y, count = _mask_and_count(seq_x, seq_y, mask, 2)

    sub_value_x = 1 + max(seq_x)
    sub_value_y = 1 + max(seq_y)

    if count == 1:
        out_x = cast(seq_x[1:])
        out_x[0] = sub_value_x
        out_y = cast(seq_y[1:])
        out_y[0] = sub_value_y
        signal = True
    else:
        out_x, out_y = cc.substitute_pairs(
            seq_x, seq_y, pair_x, pair_y, sub_value_x, sub_value_y
        )
        out_x = cast(out_x)
        out_y = cast(out_y)
        signal = False
    after = perf_counter()
    if verbose:
        return out_x, out_y, pair_x, pair_y, count, after - before, signal
    return out_x, out_y, signal


# def _onestep_windows(seq, order, verbose=True):

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


def _onestep_windows(seq, order, verbose=True):
    pass


def _onestep(seq_x, seq_y, order, verbose=True):

    if order == 2:
        return _onestep_pairs(seq_x[:], seq_y[:], verbose)
    if order > 2:
        return _onestep_windows(seq_x[:], seq_y[:], order, verbose)


def onestep(seq_x, seq_y, order, verbose=True, check=True):
    if not arraytype(seq_x):
        seq_x = cast(seq_x)
    if not arraytype(seq_y):
        seq_y = cast(seq_y)

    if check and cc.check_equality(seq_x):
        print("> All elements in sequence x are equal!")
        return None
    if check and cc.check_equality(seq_y):
        print("> All elements in sequence y are equal!")
        return None

    if len(seq_x) < order or len(seq_y) < order:
        print("> Sequence input shorter than order!\n> Can't perform substitution ...")
        return None

    return _onestep(seq_x, seq_y, order, verbose)
