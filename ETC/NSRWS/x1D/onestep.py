#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter
from itertools import compress, islice
from time import perf_counter

from ETC.NSRWS.x1D import core
from ETC.seq.recode import cast
from ETC.seq.check import arraytype


def _mask_and_count(seq, mask, order):
    """
    Apply binary mask to a sequence and count most frequently occurring windows

    This function does 3 things in the following sequence:
        1. Create sliding windows of a given size (order) - using zip and islice
        2. Apply a supplied mask to the sliding windows - using compress
        3. Count most frequently occurring window - using Counter

    In the NSRWS algorithm, this is the most time consuming step. Essentially expands
    a 1D sequence to a 2D sequence - where the sequence follows row-wise & the columnar
    expansion encodes a sliding window for each row:
        1D sequence:
            (1,2,3,4,5,6,7)

        2D expansion for window order=3:
            ((1,2,3),
             (2,3,4),
             (3,4,5),
             (4,5,6),
             (5,6,7))

        The mask is applied row-wise & must be of the same length as the number of rows
        in this 2D expansion. This is given by:
            len(mask) = len(seq) - (order - 1)

        Example application of the mask (1,0,0,1,1):
            1 -> ((1,2,3),
            0 ->  (2,3,4),    ---->      ((1,2,3),
            0 ->  (3,4,5),                (4,5,6),
            1 ->  (4,5,6),                (5,6,7))
            1 ->  (5,6,7))

        Unique windows (rows of 2D expansion) are counted and most frequently occurring
        row is returned with counts.

        1D sequence with overlap:
            (1,1,1,1,1,2,1)

        2D expansion for window order=3:
            ((1,1,1),
             (1,1,1),    ----> overlap
             (1,1,1),    ----> overlap
             (1,1,2),
             (1,2,1))

        mask will be (1,0,0,1,1) and its application will yield:
            ((1,1,1),
             (1,1,2),
             (1,2,1))

        Here, each window occurs once and the first one is returned -> (1,1,1)

    Parameters
    ----------
    seq : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    mask : array.array
        Collection of Booleans, where 0s indicate locations on "seq" to mask out.
        0s correspond to overlapping windows.
    order : int
        Size of window for NSRWS, 2 or greater.

    Returns
    -------
    freq_window : array.array
        Most frequently occurring non-overlapping "window" of size "order".
    count : int
        Number of times the most frequently occurring window occurs.

    """

    # Create overlapped sliding windows (each window a tuple of size order) & apply mask
    filtered = compress(zip(*(islice(seq, i, None) for i in range(order))), mask)

    # Count sliding windows (tuples are hashable!) & get the one most common with counts
    freq_window, count = Counter(filtered).most_common(1)[0]

    # Assign array type and return
    freq_window = cast(freq_window)

    return freq_window, count


def _onestep_pairs(seq, verbose=True):
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
    seq : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping pair substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair : array.array
            Frequent pair substituted

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
    mask = core.get_mask_pairs(seq)

    # Apply mask and find most frequent pair
    freq_pair, count = _mask_and_count(seq, mask, 2)

    # Get value for substitution of the most frequent pair with
    sub_value = 1 + max(seq)

    # If all distinct pairs, substitute the first one & set signal to True
    if count == 1:
        out = cast(seq[1:])
        out[0] = sub_value
        signal = True
    # Else, substitute all instances of the frequent pair
    else:
        out = cast(core.substitute_pairs(seq, freq_pair, sub_value))

    # Completion timer
    after = perf_counter()

    # If verbose, return more things
    if verbose:
        time_taken = after - before
        return out, signal, freq_pair, count, time_taken

    # Else return bare essentials
    return out, signal


def _onestep_windows(seq, order, verbose=True):
    """
    Execute one full step of NSRWS with order>=2 for a given sequence

    Makes use of 2 functions written in Cython & _mask_and_count in the following steps:
        1. Find overlapping windows & store their indices as mask -> get_mask_windows()
        2. Apply the mask and find most frequent window -> _mask_and_count()
        3. Substitute all occurrences of most frequent window -> substitute_windows()

    This function is different from _onestep_pairs because:
        1. This is slower due to more nested loops and checks
        2. Of course, it handles the generalized case for different window orders
        3. For higher window orders, correctness needs to be proved outside of tests

    The implementation will benefit from:
        1. Decorators for timing
        2. Decorators for verbosity of output
        3. Cython implementation of the slowest part: _mask_and_count
            problem: counting windows in C?

    Parameters
    ----------
    seq : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    order : int
        Size of window for NSRWS, 2 or greater.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping window substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair : array.array
            Frequent window substituted

        count : int
            Number of times the frequent window occurred in the sequence

        time_taken : float
            Time taken to execute step


    """

    # Initialize timer
    before = perf_counter()

    # Initialize signal for tracking sequence state with all distinct windows
    signal = False

    # Compute mask for overlapping windows
    mask = core.get_mask_windows(seq, order)

    # Apply mask and find most frequent window
    freq_window, count = _mask_and_count(seq, mask, order)

    # Get value for substitution of the most frequent window with
    sub_value = 1 + max(seq)

    # If all distinct windows, substitute the first one & set signal to True
    if count == 1:
        out = cast(seq[order - 1 :])
        out[0] = sub_value
        signal = True
    # Else, substitute all instances of the frequent window
    else:
        out = cast(core.substitute_windows(seq, order, freq_window, sub_value))

    # Completion timer
    after = perf_counter()

    # If verbose, return more things
    if verbose:
        return out, signal, freq_window, count, after - before

    # Else return bare essentials
    return out, signal


def _onestep(seq, order, verbose=True):
    """
    Wrapper that switches routine (pairs vs windows) depending on order

    For pairs (order=2), execute _onestep_pairs which is faster
    For higher orders, execute _onestep_windows

    Parameters
    ----------
    seq : array.array
        Discrete symbolic sequence containing 32-bit unsigned integers.
    order : int
        Size of window for NSRWS, 2 or greater.
    verbose : bool, optional
        Whether to report extra details. These include the frequent pair that was
        substituted, its counts & total time taken. The default is True.

    Returns
    -------
    tuple, of the following fixed elements:
        seq : array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping window substituted.

        signal : bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        freq_pair : array.array
            Frequent window substituted

        count : int
            Number of times the frequent window occurred in the sequence

        time_taken : float
            Time taken to execute step

    """

    if order == 2:
        return _onestep_pairs(seq[:], verbose)

    if order > 2:
        return _onestep_windows(seq[:], order, verbose)


def onestep(seq, order, verbose=True, check=True):
    """
    Execute one step of NSRWS on given sequence and window size.

    This function exposes the functionality of NSRWS with various checks for inputs and
    sizes. Wraps around _onestep & for convenience, allows disabling of equality check.

    Parameters
    ----------
    seq : array.array
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
    tuple, of the following fixed elements in this order:
        array.array
            Discrete symbolic sequence containing 32-bit unsigned integers, with most
            frequently occurring non-sequentially overlapping window substituted.

        bool
            indicator for the state of sequence with all distinct pairs (count=1)

    optional elements of tuple that depend on verbosity:
        array.array
            Frequent window substituted

        int
            Number of times the frequent window occurred in the sequence

        float
            Time taken to execute step

    """

    # Coerce input to appropriate array type, if not possible throw a fit & exit
    if not arraytype(seq):
        seq = cast(seq)
        if seq is None:
            return None

    # Check whether all elements are equal, if requested, & exit if True
    if check and core.check_equality(seq):
        print("> All elements in sequence are equal!")
        return None

    # Check if size of sequence is shorter than order, exit if True
    if len(seq) < order:
        print("> Sequence input shorter than order!\n> Can't perform substitution ...")
        return None

    # Else execute one step of NSRWS and return
    return _onestep(seq, order, verbose)
