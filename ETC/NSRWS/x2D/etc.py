#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from array import array

from ETC.seq import estimates as ce
from ETC.seq.recode import cast
from ETC.seq.IO import save
from ETC.NSRWS.x2D import core as cc
from ETC.NSRWS.x2D.onestep import _onestep


def _compute_verbose_truncated(seq_x, seq_y, order=2):
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

    signal = False

    # Append estimates for original sequence
    output.append(
        {
            "step": etc,
            "length": len(seq_x),
            "entropy_x": ce.entropy(seq_x),
            "entropy_y": ce.entropy(seq_y),
            "window_x": None,
            "window_y": None,
            "count": None,
            "time": None,
        }
    )

    if cc.check_equality(seq_x, seq_y):
        return etc, output

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while not signal and len(seq_x) >= order and not cc.check_equality(seq_x, seq_y):

        # Run one step of NSRWS in verbose mode (returns window and count)
        seq_x, seq_y, signal, pair_x, pair_y, count, time = _onestep(
            seq_x, seq_y, order, verbose=True
        )

        # Increment ETC
        etc += 1

        # Compute estimates and append to aggregator
        output.append(
            {
                "step": etc,
                "length": len(seq_x),
                "entropy_x": ce.entropy(seq_x),
                "entropy_y": ce.entropy(seq_y),
                "window_x": pair_x,
                "window_y": pair_y,
                "count": count,
                "time": time,
            }
        )
    n = 0
    if signal and not cc.check_equality(seq_x, seq_y):

        while len(seq_x) >= order and n < 5:
            # Run one step of NSRWS in verbose mode (returns window and count)
            seq_x, seq_y, signal, pair_x, pair_y, count, time = _onestep(
                seq_x, seq_y, order, verbose=True
            )

            # Increment ETC
            etc += 1

            # Compute estimates and append to aggregator
            output.append(
                {
                    "step": etc,
                    "length": len(seq_x),
                    "entropy_x": ce.entropy(seq_x),
                    "entropy_y": ce.entropy(seq_y),
                    "window_x": pair_x,
                    "window_y": pair_y,
                    "count": count,
                    "time": time,
                }
            )
            n += 1
        if len(seq_x) % (order - 1) == 0:
            etc += len(seq_x) // (order - 1) - 1
        else:
            etc += len(seq_x) // (order - 1)

    # Display ETC and return it with aggregator
    # print(f"ETC={etc}")
    return etc, output


def _compute_verbose_full(seq_x, seq_y, order=2):
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
            "length": len(seq_x),
            "entropy_x": ce.entropy(seq_x),
            "entropy_y": ce.entropy(seq_y),
            "window_x": None,
            "window_y": None,
            "count": None,
            "time": None,
        }
    )

    if cc.check_equality(seq_x, seq_y):
        return etc, output

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while len(seq_x) >= order and not cc.check_equality(seq_x, seq_y):

        # Run one step of NSRWS in verbose mode (returns window and count)
        seq_x, seq_y, signal, pair_x, pair_y, count, time = _onestep(
            seq_x, seq_y, order, verbose=True
        )
        # Increment ETC
        etc += 1

        # Compute estimates and append to aggregator
        output.append(
            {
                "step": etc,
                "length": len(seq_x),
                "entropy_x": ce.entropy(seq_x),
                "entropy_y": ce.entropy(seq_y),
                "window_x": pair_x,
                "window_y": pair_y,
                "count": count,
                "time": time,
            }
        )
    # Display ETC and return it with aggregator
    # print(f"ETC={etc}")
    return etc, output


def _compute_compact_truncated(seq_x, seq_y, order=2):
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

    if cc.check_equality(seq_x, seq_y):
        return etc

    signal = False
    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    # while len(seq) >= order and not equality(seq):
    while not signal and len(seq_x) >= order and not cc.check_equality(seq_x, seq_y):

        # Run one step of NSRWS
        seq_x, seq_y, signal = _onestep(seq_x, seq_y, order, verbose=False)

        # Increment ETC
        etc += 1

    n = 0
    if signal and not cc.check_equality(seq_x, seq_y):

        while len(seq_x) >= order and n < 5:
            # Run one step of NSRWS in verbose mode (returns window and count)
            seq_x, seq_y, signal = _onestep(seq_x, seq_y, order, verbose=False)

            # Increment ETC
            etc += 1

            n += 1

        if len(seq_x) % (order - 1) == 0:
            etc += len(seq_x) // (order - 1) - 1
        else:
            etc += len(seq_x) // (order - 1)

    # Display ETC and return it
    # print(f"ETC={etc}")
    return etc


def _compute_compact_full(seq_x, seq_y, order=2):
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

    if cc.check_equality(seq_x, seq_y):
        return etc

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    # while len(seq) >= order and not equality(seq):
    while len(seq_x) >= order and not cc.check_equality(seq_x, seq_y):

        # Run one step of NSRWS
        seq_x, seq_y, signal = _onestep(seq_x, seq_y, order, verbose=False)

        # Increment ETC
        etc += 1

    # Display ETC and return it
    # print(f"ETC={etc}")
    return etc


def compute(seq_x, seq_y, order=2, verbose=True, truncate=True):
    """
    This function estimates the Effort-To-Compress for a given sequence. It
    wraps around other functions and executes them based on input options.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    order : int, optional
        Number of elements in window for substitution.
        The default is 2 for pairs.
    verbose : bool, optional
        Whether to compute additional metrics. The default is True.

    Returns
    -------
    dict
        ETC (int) & optionally, trajectory of NSRWS algorithm if verbose=True.

    """
    assert len(seq_x) == len(seq_y), "ERROR: The 2 sequences should have the same length!"

    # Create a copy of the original sequence
    seq_x = cast(seq_x)
    seq_y = cast(seq_y)

    if truncate:
        # If verbose, run the verbose version and return accordingly
        if verbose:
            etc, out = _compute_verbose_truncated(seq_x, seq_y, order)
            return {"ETC2D": etc, "NETC2D": etc / (len(seq_x) - 1), "Trajectory": out}
        else:
            # If not verbose, run the compact version and return accordingly
            etc = _compute_compact_truncated(seq_x, seq_y, order)
            return {"ETC2D": etc, "NETC2D": etc / (len(seq_x) - 1)}
    else:
        # If verbose, run the verbose version and return accordingly
        if verbose:
            etc, out = _compute_verbose_full(seq_x, seq_y, order)
            return {"ETC2D": etc, "NETC2D": etc / (len(seq_x) - 1), "Trajectory": out}
        else:
            # If not verbose, run the compact version and return accordingly
            etc = _compute_compact_full(seq_x, seq_y, order)
            return {"ETC2D": etc, "NETC2D": etc / (len(seq_x) - 1)}


def compute_save(seq_x, seq_y, filename, truncate=True, order=2):
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

    Returns
    -------
    dict
        ETC (int).

    """
    assert len(seq_x) == len(seq_y), "ERROR: The 2 sequences should have the same length!"

    # Create a copy of the original sequence
    seq_x = cast(seq_x)
    seq_y = cast(seq_y)

    if truncate:
        etc, out = _compute_verbose_truncated(seq_x, seq_y, order)
    else:
        etc, out = _compute_verbose_full(seq_x, seq_y, order)

    # Save the output to a csv file and return
    save(out, filename)

    return {"ETC2D": etc, "NETC2D": etc / (len(seq_x) - 1)}


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
