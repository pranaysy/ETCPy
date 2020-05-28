#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains a generalized implementation of the NSRPS algorithm, from
here on referred to as the NSRWS: Non-Sequential Recursive Window Subsitution.

Given an input sequence of integer values, as a list/tuple, it returns this
sequence with the most frequently occurring non-sequential window substituted.

The module consists of 5 modular sub-functions that carry out multiple steps
of the algorithm. The function run_once_NSRWS wraps around all sub-functions &
executes the full algorithm once.

External Dependecies: None

@author: Pranay S. Yadav
"""

from collections import Counter

# Import functions from standard library modules
from itertools import compress, islice, repeat

# Import functions from local modules
from ETC.utils import equality


# Function definitions
def _find_overlapping_windows(seq, order=2):
    """
    This function takes in a collection (list or tuple) of integers and returns
    a tuple of overlapping sliding windows (shifted by one element) each of size
    'order' elements. Total number of windows returned depends on order and is
    (sequence length - order + 1)

    Makes use of islice from itertools, and is memory efficient by avoiding
    creating copies of the input sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    windows : tuple
        Tuple of tuples, each element is a sliding window of size=order.

    """
    # Use islice for generating sliding cascades from the sequence, iteratively
    windows = tuple(zip(*(islice(seq, i, None) for i in range(order))))

    return windows


def _filter_pairs(pairs):
    """
    This function takes a tuple of 2-tuples (output of _find_overlapping_windows,
    with order=2) and returns a list of booleans representing a mask for windows.
    This mask has False elements at indices corresponding to sequentially
    overlapping pairs which are identical.

    Example:
        for an input sequence: (2,1,1,1,1,2)
        win = ((2,1), (1,1), (1,1), (1,1), (1,2))
        mask = [True, True, False, True, True]
        * 3rd window overlaps with the 2nd and is masked out, but 4th is kept

    Parameters
    ----------
    pairs : tuple
        Tuple of 2-tuples, each element is sliding pair of elements.

    Returns
    -------
    mask : list
        Selection mask for pairs where False elements correspond to those
        successive overlapping pairs which are identical.

    """
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


def _filter_windows(win, order):
    """
    This function takes a tuple of tuples (output of _find_overlapping_windows)
    and returns a list of booleans representing a mask for windows. This mask
    has False elements at indices corresponding to sequentially overlapping
    windows which are identical. Sequential non-overlapping windows that are
    identical are preserved.

    The number of sequential overlapping windows is one less than 'order'. In
    other words, for order=n and for sequentially identical windows, every nth
    window will be selected beginning from the first of these windows.

    Example A for order=3:
        for an input sequence: (2,1,1,1,1,2)
        win = ((2,1,1), (1,1,1), (1,1,1), (1,1,2))
        mask = [True, True, False, True]
        * 3rd window overlaps with the 2nd and is masked out

    Example B for order=3:
        for an input sequence: (2,1,1,1,1,1,1,1,2)
        win = ((2,1,1), (1,1,1), (1,1,1), (1,1,1), (1,1,1), (1,1,1), (1,1,2))
        mask = [True, True, False, False, True, False, True]
        * 3rd, 4th, 5th & 6th windows are identical to 2nd, but 3rd & 4th
        * share overlap with 2nd & are filtered out, but 5th is not. Similarly,
        * 6th overlaps with 5th and is filtered out.

    Parameters
    ----------
    win : tuple
        Tuple of tuples, each element is a sliding window of size=order.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    mask : list
        Selection mask for windows where False elements correspond to those
        successive overlapping windows which are identical.

    """
    # Initialize mask to all True
    mask = list(repeat(True, len(win)))

    # Loop begins at 0, ends on penultimate index (each check looks ahead by 1)
    idx = 0
    upper_idx_lim = len(win) - order + 1

    # Main while loop
    while idx < upper_idx_lim:

        # Check equality for all remaining overlapping windows
        ## For order 'n', there are 'n-1' sequentially overlapping windows
        ## As the 1st pair has already been compared, 'n-2' windows remain
        for subidx in range(idx, idx + order - 1):
            counter = False
            print(idx, subidx + 1, win[idx], win[subidx + 1])
            # Check if index is valid, and if so, check if window at
            # current & next index are equal
            if win[idx] == win[subidx + 1]:

                # If so, next window is masked out, current window is "kept in"
                # mask[subidx + 1] = False

                # Increment the overall index to slide by 1
                counter = True

        if counter:
            for n in range(order - 1):
                mask[idx + 1 + n] = False
        print("\t", mask, "\n\t", win)
        # Increment while loop index
        idx += order

    # Return the mask
    return mask


def _apply_filter_mask(win, mask):
    """
    This function applies a boolean mask for selection of windows from a tuple
    of windows. This is used for selecting sequentially non-overlapping windows.
    It also applies a mask to the indices for preserving order and locations.

    Makes use of compress from itertools and returns 2 tuples, one for masked
    tuple of windows and another for masked indices.

    Parameters
    ----------
    win : tuple
        Tuple of tuples, each element is a sliding window.
    mask : list
        Selection mask for windows where False elements correspond to those
        successive overlapping windows which are identical. Should have same
        length as win.

    Returns
    -------
    filtered_windows : tuple
        Tuple of tuples, each element is a sliding window with sequentially
        identical overlapping windows removed.
    filtered_indices : tuple
        Indices corresponding to original locations of filtered_windows in win.

    """
    # Apply mask to win
    filtered_windows = tuple(compress(win, mask))

    # Generate indices and apply mask to them
    filtered_indices = tuple(compress(range(len(win)), mask))

    return filtered_windows, filtered_indices


def _find_frequent_windows(filtered_windows, filtered_indices):
    """
    This function returns the most frequently occurring window from a tuple of
    windows along with its counts and indices in the tuple.

    Relies on behavior of the most_common() method of Counter objects, where if
    multiple windows with same counts are found, then the one first encountered
    is returned.

    Parameters
    ----------
    filtered_windows : tuple
        Tuple of tuples, each element is a sliding window with sequentially
        identical overlapping windows removed.
    filtered_indices : tuple
        Indices corresponding to original locations of filtered_windows in win.

    Returns
    -------
    freq_win : tuple
        Most frequently occurring window. If multiple windows with equal counts,
        then this is the first window encountered amongst them.
    count : int
        Number of times the most frequent window occurred.
    idx_freq_win : compress iterator
        Indices that correspond to the occurence of the most frequent window in
        filtered_windows. Since the filtered_indices correspond to the original
        locations of windows, these indices directly correspond to them too.

    """
    # Count all unique windows, return the *one* most common window
    freq_win, count = Counter(filtered_windows).most_common(1)[0]

    # Select only those indices that correspond to the most frequent window
    idx_freq_win = compress(filtered_indices, [a == freq_win for a in filtered_windows])

    return freq_win, count, idx_freq_win


def _substitute_window(seq, idx_freq_win, order=2):
    """
    This function carries out substitution of the most frequent window with an
    integer. Each window substitution reduces the length of the sequence by
    (order-1) elements. The substituted integer is one greater than the largest
    element in the sequence.

    The substitution uses indices/locations of the most frequent window and
    is carried out in 3 steps, for a single window:
        1. Replace 1st element of seq at window index by integer (max(seq)+1)
        2. Replace the next order - 1 elements with False
        3. Filter list to eliminate all False elements
    Iteration of first 2 steps happens together over indices in idx_freq_win

    Parameters
    ----------
    seq : list
        Sequence of integers.
    idx_freq_win : compress iterator or tuple or list
        Indices that correspond to the occurence of the most frequent window.
        The index of the window is the same as the index of the first element
        of that window in seq.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    reduced_seq : list
        Sequence of integers with most frequent window substituted.
            len(reduced_seq) < len(seq)

    """
    # Integer value to use for substitution
    a = max(seq) + 1

    # Iterate over indices/locations of most frequent window
    for index in idx_freq_win:

        # Substitute the first element of window, in seq, by integer
        seq[index] = a

        # Substitute the remaining elements of window, in seq, with False
        for n in range(1, order):
            seq[index + n] = False

    # Filter out all False elements and return what's left
    reduced_seq = [element for element in seq if element]

    return reduced_seq


def _execute_onestep(seq, order=2, verbose=False):
    """
    This function runs one full step of the non-sequential recursive window-
    substitution algorithm. (For pairs or order=2, called NSRPS).

    For internal use only, as this function does not carry out sanity checks.
    For general/external usage, refer to run_once_NSRPS.

    Parameters
    ----------
    seq : list
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.
    verbose : bool, optional
        If True, returns most frequent pair with counts. The default is False.

    Returns
    -------
    reduced_seq : list
        Sequence of integers with most frequent window substituted.
    The following are returned if verbose=True
        freq_win : tuple
            Most frequent window in seq.
        count : int
            Number of times the most frequent window occurred in seq.

    """
    # Get sliding overlapping windows
    windows = _find_overlapping_windows(seq, order)

    if order == 2:
        # Get mask that filters overlapping pairs out
        mask = _filter_pairs(windows)
    else:
        # Get mask that filters overlapping windows out
        mask = _filter_windows(windows, order)

    # Apply mask, get filtered windows and indices
    filt_win, filt_idx = _apply_filter_mask(windows, mask)

    # Get most frequent window, with counts and indices
    freq_win, count, idx_freq_win = _find_frequent_windows(filt_win, filt_idx)

    # Carry out substitution of most frequent window at given indices
    reduced_seq = _substitute_window(seq, idx_freq_win, order)

    # Return substituted sequence with optional outputs
    if verbose:
        return reduced_seq, freq_win, count
    return reduced_seq


def run_once_NSRWS(seq, order=2, check=True, verbose=False):
    """
    This function executes one step of the Non-Sequential Recursive Window
    Substitution algorithm. Defaults to a window size (order) of 2, called
    Non-Sequential Recursive Pair Substitution or NSRPS.

    Uses equality from ETC.utils

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.
    check : bool, optional
        If True, checks if all elements are same. The default is True.
    verbose : bool, optional
        If True, returns most frequent pair with counts. The default is False.

    Returns
    -------
    reduced_seq : list
        Sequence of integers with most frequent window substituted.
    The following are returned if verbose=True
        freq_win : tuple
            Most frequent window in seq.
        count : int
            Number of times the most frequent window occurred in seq.

    """
    # Convert collection to mutable list & work on a copy to avoid side-effects
    if isinstance(seq, tuple):
        temp = list(seq)
    else:
        temp = seq.copy()

    # If check is True, return nothing if all elements are equal
    if check:
        if equality(seq):
            print("> Constant symbolic sequence, nothing to substitute ...")
            return None
        else:
            return _execute_onestep(temp, order, verbose)

    # If check is disabled, run one step and return
    return _execute_onestep(temp, order, verbose)
