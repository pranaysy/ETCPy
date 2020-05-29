#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter
from itertools import compress, islice
from time import perf_counter

from ETC.NSRWS.x2D import core
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
    filtered = compress(
        zip(*(islice(zip(seq_x, seq_y), i, None) for i in range(order))), mask
    )

    # Count sliding windows (tuples are hashable!) & get the one most common with counts
    freq_pair, count = Counter(filtered).most_common(1)[0]

    # Assign array type and return
    pair_x = cast([freq_pair[0][0], freq_pair[1][0]])
    pair_y = cast([freq_pair[0][1], freq_pair[1][1]])

    return pair_x, pair_y, count


def _onestep_pairs(seq_x, seq_y, verbose=True):
    """
    Execute one full step of NSRPS (NSRWS with order=2) for a given sequence

    Makes use of 2 functions written in Cython & _mask_and_count in the following steps:
        1. Find overlapping pairs & store their indices for masking -> get_mask_pairs()
        2. Apply the mask and find most frequent pair -> _mask_and_count()
        3. Substitute all occurrences of the most frequent pair -> substitute_pairs()

    This function is different from _onestep_windows because:
        1. It is *much* faster due to fewer nested loops
        2. It targets a more common use case scenario: for distances, for CCC, etc
        3. For higher window orders, correctness needs to be proved outside of tests

    The implementation will benefit from:
        1. Decorators for timing
        2. Decorators for verbosity of output
        3. Cython implementation of the slowest part: _mask_and_count
            problem: counting windows in C?

    Parameters
    ----------
    seq_x : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    seq_y : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq_x : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        seq_y : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair_x : array.array
            Frequent pair substituted in seq_x

        freq_pair_y : array.array
            Frequent pair substituted in seq_y

        count : int
            Number of times the frequent pair occurred in the sequence

        time_taken : float
            Time taken to execute step


    """
    # Initialize timer
    before = perf_counter()

    # Initialize signal for tracking sequence state with all distinct pairs
    signal = False

    # Compute mask for overlapping pairs
    mask = core.get_mask_pairs(seq_x, seq_y)

    # Apply mask and find most frequent pair
    pair_x, pair_y, count = _mask_and_count(seq_x, seq_y, mask, 2)

    # Get values for substitution of the most frequent pair with
    sub_value_x = 1 + max(seq_x)
    sub_value_y = 1 + max(seq_y)

    # If all distinct pairs, substitute the first one & set signal to True
    if count == 1:
        out_x = cast(seq_x[1:])
        out_x[0] = sub_value_x

        out_y = cast(seq_y[1:])
        out_y[0] = sub_value_y

        signal = True
    # Else, substitute all instances of the frequent pair
    else:
        out_x, out_y = core.substitute_pairs(
            seq_x, seq_y, pair_x, pair_y, sub_value_x, sub_value_y
        )
        out_x = cast(out_x)
        out_y = cast(out_y)

    # Completion timer
    after = perf_counter()

    # If verbose, return more things
    if verbose:
        time_taken = after - before
        return out_x, out_y, signal, pair_x, pair_y, count, time_taken

    # Else return bare essentials
    return out_x, out_y, signal


# def _onestep_windows(seq, order, verbose=True):

#     before = perf_counter()
#     mask = core.get_mask_windows(seq, order)[: -(order - 1)]
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
#         out = array("I", core.substitute_windows(seq, order, window, sub_value))
#         signal = False
#     out = array("I", core.substitute_windows(seq, order, window, sub_value))
#     after = perf_counter()
#     if verbose:
#         return out, freq_window, count, after - before, signal
#     return out, signal


# def _onestep_windows(seq, order, verbose=True):
#     pass


def _onestep(seq_x, seq_y, order, verbose=True):
    """
    Wrapper that switches routine (pairs vs windows) depending on order

    For pairs (order=2), execute _onestep_pairs which is faster
    For higher orders, execute _onestep_windows

    Parameters
    ----------
    seq_x : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    seq_y : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    order : int
        Size of window for NSRWS, 2 or greater.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq_x : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        seq_y : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair_x : array.array
            Frequent pair substituted in seq_x

        freq_pair_y : array.array
            Frequent pair substituted in seq_y

        count : int
            Number of times the frequent pair occurred in the sequence

        time_taken : float
            Time taken to execute step

    """
    if order == 2:
        return _onestep_pairs(seq_x[:], seq_y[:], verbose)
    # if order > 2:
    #     return _onestep_windows(seq_x[:], seq_y[:], order, verbose)


def onestep(seq_x, seq_y, order, verbose=True, check=True):
    """
    Execute one step of NSRWS on given sequence and window size.

    This function exposes the functionality of NSRWS with various checks for inputs and
    sizes. Wraps around _onestep & for convenience, allows disabling of equality check.

    Parameters
    ----------
    seq_x : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    seq_y : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    order : int
        Size of window for NSRWS, 2 or greater.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.
    check : bool, optional
        Check for equality of all symbols in sequence. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq_x : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        seq_y : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair_x : array.array
            Frequent pair substituted in seq_x

        freq_pair_y : array.array
            Frequent pair substituted in seq_y

        count : int
            Number of times the frequent pair occurred in the sequence

        time_taken : float
            Time taken to execute step

    """
    # Check if both sequences are of same length, if not then exit
    if len(seq_x) != len(seq_y):
        print("> Both inputs must be of the same length!")
        return None

    # Coerce input 1 to appropriate array type, if not possible throw a fit & exit
    if not arraytype(seq_x):
        seq_x = cast(seq_x)

    # Coerce input 2 to appropriate array type, if not possible throw a fit & exit
    if not arraytype(seq_y):
        seq_y = cast(seq_y)

    # Exit if neither inputs could be coerced
    if seq_x is None or seq_y is None:
        return None

    # Check if size of sequence is shorter than order, exit if True
    if len(seq_x) < order or len(seq_y) < order:
        print("> Sequence input shorter than order!\n> Can't perform substitution ...")
        return None

    # Check whether all elements are equal, if requested, & exit if True
    if check and core.check_equality(seq_x, seq_y):
        print("> All elements in sequence x are equal!")
        return None

    # Else execute one step of NSRWS and return
    return _onestep(seq_x, seq_y, order, verbose)
