#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from ETC.seq import estimates as ce
from ETC.seq.recode import cast
from ETC.seq.IO import save
from ETC.NSRWS.x1D import core as cc
from ETC.NSRWS.x1D.onestep import _onestep


def _compute_verbose_truncated(seq, order=2):
    """
    This function runs the NSRWS algorithm for estimation of ETC and extracts
    additional metrics at each step of the algorithm. These include:
        - length of sequence
        - entropy of sequence
        - most frequent window
        - count of most frequent window

    The NSRWS algorithm is run iteratively until all elements are equal or the
    sequence has been reduced to a size smaller than the size of the window
    being substituted (specified by order). The number of steps taken till the
    iteration stops is the Effort-To-Compress (ETC) estimate for the sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    etc : int
        Effort-To-Compress estimate for given seq and order.
    output : list
        List of dictionaries corresponding to each step of NSRWS run during
        estimation of ETC for the given sequence.

    """
    # Initialize ETC to 0
    etc = 0

    # Initialize an aggregator for collecting dictionaries of estimates
    output = list()

    # Append estimates for original sequence
    output.append(
        {
            "step": etc,
            "length": len(seq),
            "entropy": ce.entropy(seq),
            "window": None,
            "count": None,
            "time": None,
        }
    )

    # Check if all elements are equal and break if so
    if cc.check_equality(seq):
        return etc, output

    # Initialize a boolean for tracking truncation step
    signal = False

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while not signal and len(seq) >= order and not cc.check_equality(seq):

        # Run one step of NSRWS in verbose mode (returns window and count)
        seq, signal, freq_win, count, time = _onestep(seq, order, verbose=True)

        # Increment ETC
        etc += 1

        # Compute estimates and append to aggregator
        output.append(
            {
                "step": etc,
                "length": len(seq),
                "entropy": ce.entropy(seq),
                "window": freq_win,
                "count": count,
                "time": time,
            }
        )

    if signal:  # Run 5 times for estimating entropy
        n = 0
        while len(seq) >= order and not cc.check_equality(seq) and n < 5:

            # Run one step of NSRWS in verbose mode (returns window and count)
            seq, signal, freq_win, count, time = _onestep(seq, order, verbose=True)

            # Increment ETC
            etc += 1

            # Compute estimates and append to aggregator
            output.append(
                {
                    "step": etc,
                    "length": len(seq),
                    "entropy": ce.entropy(seq),
                    "window": freq_win,
                    "count": count,
                    "time": time,
                }
            )

            # Increment while loop indexer
            n += 1

        # Truncation logic
        if len(seq) % (order - 1) == 0:  # If no remainder
            etc += len(seq) // (order - 1) - 1  # one less step is needed
        else:
            etc += len(seq) // (order - 1)
    # Return ETC it with aggregator
    return etc, output


def _compute_verbose_full(seq, order=2):
    """
    This function runs the NSRWS algorithm for estimation of ETC and extracts
    additional metrics at each step of the algorithm. These include:
        - length of sequence
        - entropy of sequence
        - most frequent window
        - count of most frequent window

    The NSRWS algorithm is run iteratively until all elements are equal or the
    sequence has been reduced to a size smaller than the size of the window
    being substituted (specified by order). The number of steps taken till the
    iteration stops is the Effort-To-Compress (ETC) estimate for the sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    etc : int
        Effort-To-Compress estimate for given seq and order.
    output : list
        List of dictionaries corresponding to each step of NSRWS run during
        estimation of ETC for the given sequence.

    """
    # Initialize ETC to 0
    etc = 0

    # Initialize an aggregator for collecting dictionaries of estimates
    output = list()

    # Append estimates for original sequence
    output.append(
        {
            "step": etc,
            "length": len(seq),
            "entropy": ce.entropy(seq),
            "window": None,
            "count": None,
            "time": None,
        }
    )

    # Check if all elements are equal and break if so
    if cc.check_equality(seq):
        return etc, output

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while len(seq) >= order and not cc.check_equality(seq):

        # Run one step of NSRWS in verbose mode (returns window and count)
        seq, signal, freq_win, count, time = _onestep(seq, order, verbose=True)

        # Increment ETC
        etc += 1

        # Compute estimates and append to aggregator
        output.append(
            {
                "step": etc,
                "length": len(seq),
                "entropy": ce.entropy(seq),
                "window": freq_win,
                "count": count,
                "time": time,
            }
        )
    # Return ETC with aggregator
    return etc, output


def _compute_compact_truncated(seq, order=2):
    """
    This function runs the NSRWS algorithm for estimation of ETC. This is a
    compact version of the _compute_verbose function and does not return any
    additional metrics.

    The NSRWS algorithm is run iteratively until all elements are equal or the
    sequence has been reduced to a size smaller than the size of the window
    being substituted (specified by order). The number of steps taken till the
    iteration stops is the Effort-To-Compress (ETC) estimate for the sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    etc : int
        Effort-To-Compress estimate for given seq and order.

    """
    # Initialize ETC to 0
    etc = 0

    # Check if all elements are equal and break if so
    if cc.check_equality(seq):
        return etc

    # Initialize a boolean for tracking truncation step
    signal = False

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while not signal and len(seq) >= order and not cc.check_equality(seq):

        # Run one step of NSRWS
        seq, signal = _onestep(seq, order, verbose=False)

        # Increment ETC
        etc += 1

    if signal:

        # Run NSRWS twice, windows will definitely become unique in 2 steps?
        n = 0
        while len(seq) >= order and n < 2:

            # Run one step of NSRWS in verbose mode (returns window and count)
            seq, signal = _onestep(seq, order, verbose=False)

            # Increment ETC
            etc += 1

            # Increment while loop indexer
            n += 1

        # Truncation logic
        if len(seq) % (order - 1) == 0:  # If no remainder
            etc += len(seq) // (order - 1) - 1  # One less step is needed
        else:
            etc += len(seq) // (order - 1)

    # Display ETC and return it
    # print(f"ETC={etc}")
    return etc


def _compute_compact_full(seq, order=2):
    """
    This function runs the NSRWS algorithm for estimation of ETC. This is a
    compact version of the _compute_verbose function and does not return any
    additional metrics.

    The NSRWS algorithm is run iteratively until all elements are equal or the
    sequence has been reduced to a size smaller than the size of the window
    being substituted (specified by order). The number of steps taken till the
    iteration stops is the Effort-To-Compress (ETC) estimate for the sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.

    Returns
    -------
    etc : int
        Effort-To-Compress estimate for given seq and order.

    """
    # Initialize ETC to 0
    etc = 0

    # Check if all elements are equal and break if so
    if cc.check_equality(seq):
        return etc

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while len(seq) >= order and not cc.check_equality(seq):

        # Run one step of NSRWS
        seq, signal = _onestep(seq, order, verbose=False)

        # Increment ETC
        etc += 1

    # Display ETC and return it
    # print(f"ETC={etc}")
    return etc


def compute(seq, order=2, verbose=False, truncate=True):
    """
    Estimate the Effort-To-Compress for a given sequence using the NSRPS algorithm.

    This function wraps around other functions and switches between them based on input
    parameters. The default options give the fastest results.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.
    verbose : bool, optional
        Whether to compute additional metrics. The default is True.
    truncate: bool, optional
        Whether to halt iterative estimation once fully saturated 'axiom' has been reached

    Returns
    -------
    dict
        ETC1D (int), NETC1D (float) & optionally, trajectory of algorithm if verbose=True

    """
    # Create a copy of the original sequence with the appropriate type
    seq = cast(seq)

    if truncate:
        # If verbose, run the verbose version and return accordingly
        if verbose:
            etc, out = _compute_verbose_truncated(seq, order)
            return {"ETC1D": etc, "NETC1D": etc / (len(seq) - 1), "Trajectory": out}
        else:
            # If not verbose, run the compact version and return accordingly
            etc = _compute_compact_truncated(seq, order)
            return {"ETC1D": etc, "NETC1D": etc / (len(seq) - 1)}
    else:
        # If verbose, run the verbose version and return accordingly
        if verbose:
            etc, out = _compute_verbose_full(seq, order)
            return {"ETC1D": etc, "NETC1D": etc / (len(seq) - 1), "Trajectory": out}
        else:
            # If not verbose, run the compact version and return accordingly
            etc = _compute_compact_full(seq, order)
            return {"ETC1D": etc, "NETC1D": etc / (len(seq) - 1)}


def compute_save(seq, filename, order=2, truncate=True):
    """
    This function estimates the Effort-To-Compress for a given sequence in
    verbose mode and writes the trajectory of the NSRWS algorithm to disk.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    filename : str or Path object
        Name of output file or path to output file.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.
    truncate: bool, optional
        Whether to halt iterative estimation once fully saturated 'axiom' has been reached

    Returns
    -------
    dict
        ETC1D (int), NETC1D (float).

    """
    # Create a copy of the original sequence with the appropriate type
    seq = cast(seq)

    if truncate:
        etc, out = _compute_verbose_truncated(seq, order)
    else:
        etc, out = _compute_verbose_full(seq, order)

    # Save the output to a csv file and return
    save(out, filename)

    return {"ETC1D": etc, "NETC1D": etc / (len(seq) - 1)}


# %% Adapt into test!
# from random import randint
# iterations = 100
# misses = 0
# for _ in range(iterations):
#     n = randint(50,2000)
#     q = range(n)
#     order = 2
#     if n%(order-1) == 0:
#         etc = n//(order-1) -1
#     else:
#         etc = n//(order-1)
#     out = compute(q, order, 0, 0)['ETC']
#     if etc != out:
#         misses += 1
#     print(order, n, etc, out, etc==out, misses)
# print('total misses:', 100*misses/iterations)

# %% Adapt into another test!
# misses = 0
# for order in range(2, 50):
#     if compute(x, order, 0, 1)['ETC'] != compute(x, order, 0, 0)['ETC']:
#         misses+=1
#         print(order)
