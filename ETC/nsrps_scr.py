#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter, defaultdict
from itertools import compress, islice, repeat

import numpy as np

from ETC.nsrps_c import runner_count as rc
from ETC.nsrps_c import runX


def _get_windows(seq, order):
    return tuple(zip(*(islice(seq, i, None) for i in range(order))))


def _filter_pairs(pairs):
    # Initialize mask to all True
    mask = list(repeat(True, len(pairs)))

    # Loop begins at 0, ends on penultimate index (each check looks ahead by 1)
    idx = 0
    upper_idx_lim = len(pairs) - 1

    # Main while loop
    while idx < upper_idx_lim:

        # Check if current & next windows equal (tuple comparison, short-circuit)
        if pairs[idx] == pairs[idx + 1]:

            # If so, next window is masked out, current window is "kept in"
            mask[idx + 1] = False

            # Move on to the next window
            idx += 1

        # Increment while loop index
        idx += 1

    # Return the mask
    return mask


def _apply_filter_mask(win, mask):

    # Apply mask to win
    filtered_windows = tuple(compress(win, mask))

    # Generate indices and apply mask to them
    filtered_indices = tuple(compress(range(len(win)), mask))

    return filtered_windows, filtered_indices


def _find_frequent_windows(filtered_windows, filtered_indices):

    # Count all unique windows, return the *one* most common window
    freq_win, count = Counter(filtered_windows).most_common(1)[0]

    # Select only those indices that correspond to the most frequent window
    idx_freq_win = tuple(
        compress(filtered_indices, [a == freq_win for a in filtered_windows])
    )

    return freq_win, count, idx_freq_win


def old_skool(seq, order=2):
    win = _get_windows(seq, order)
    mask = _filter_pairs(win)
    filtwin, filtidx = _apply_filter_mask(win, mask)
    freq_win, count, idx_freq_win = _find_frequent_windows(filtwin, filtidx)
    return freq_win, count, idx_freq_win


def runner(x, order=2):
    # iterate over x
    n = 0
    mask = [1] * len(b)
    while n < len(x) - order + 1:

        if x[n : n + 2] == x[n + 1 : n + 3]:
            mask[n + 1] = 0
            n += 1

        n += 1

    return mask


def runner_counts(x, order=2):
    # iterate over x
    n = 0
    counts = Counter()
    idxes = defaultdict(list)

    while n < len(x) - order + 1:

        if x[n : n + order] == x[n + 1 : n + order + 1]:
            n += 1
        idxes[tuple(x[n : n + order])].append(n)
        counts[tuple(x[n : n + order])] += 1
        n += 1

    freq_win, counts = counts.most_common(1)[0]
    idx_freq_win = idxes[freq_win]

    return freq_win, counts, idx_freq_win


b = [1] * 100000
bn = np.array(b, dtype="u4")


def test(x):
    # iterate over x

    order = 2
    # cdef long int n = 0
    # cdef long int m = 0
    # cdef long int matches = 0
    upperlim = x.shape[0] - order
    n = m = 0
    # counts = Counter()
    # idxes = defaultdict(list)

    while n < upperlim:
        matches = 0
        for m in range(order):
            if x[n + m] == x[n + m + 1]:
                matches += 1
        if matches == order:
            n += 1

        # idxes[tuple(x[n:n+order])].append(n)
        # counts[tuple(x[n:n+order])] += 1
        n += 1

    # freq_win, counts = counts.most_common(1)[0]
    # idx_freq_win = idxes[freq_win]

    return n
