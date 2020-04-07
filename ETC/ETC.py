#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains an implementation of the Effort-To-Compress (ETC) algorithm
that makes iterative use of Non-Sequential Recursive Window (Pair) Substitution
as described in Nagaraj et al, 2013

The module consists of 2 main functions:
    compute - estimates ETC for a given sequence with desired verbosity
    compute_save - estimates ETC verbosely and saves output to disk

External Dependecies: None

Reference paper:
    Nagaraj, Nithin, Karthi Balasubramanian, and Sutirth Dey. A New Complexity
    Measure for Time Series Analysis and Classification. The European Physical
    Journal Special Topics 222, no. 3–4 (July 2013): 847–60.
    https://doi.org/10.1140/epjst/e2013-01888-9.

@author: Pranay S. Yadav
"""

# Import functions from local modules
from ETC.NSRWS import _execute_one_step
from ETC.utils import equality, entropy
from ETC.IO import save

# Function definitions
def _compute_verbose(seq, order=2):
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

    # Create a copy of the original sequence
    temp = list(seq)

    # Initialize an aggregator for collecting dictionaries of estimates
    output = list()

    # Append estimates for original sequence
    output.append(
        {
            "step": etc,
            "length": len(temp),
            "entropy": entropy(temp),
            "window": None,
            "count": None,
        }
    )

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while len(temp) >= order and not equality(temp):

        # Run one step of NSRWS in verbose mode (returns window and count)
        temp, freq_win, count = _execute_one_step(temp, order, verbose=True)

        # Increment ETC
        etc += 1

        # Compute estimates and append to aggregator
        output.append(
            {
                "step": etc,
                "length": len(temp),
                "entropy": entropy(temp),
                "window": freq_win,
                "count": count,
            }
        )

    # Display ETC and return it with aggregator
    print(f"ETC={etc}")
    return etc, output


def _compute_compact(seq, order=2):
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

    # Create a copy of the original sequence
    temp = list(seq)

    # Execute iteration loop until either all elements are equal or sequence is
    # reduced to less than size of the window being substituted (order)
    while len(temp) >= order and not equality(temp):

        # Run one step of NSRWS
        temp = _execute_one_step(temp, order, verbose=False)

        # Increment ETC
        etc += 1

    # Display ETC and return it
    print(f"ETC={etc}")
    return etc


def compute(seq, order=2, verbose=True):
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
    # If verbose, run the verbose version and return accordingly
    if verbose:
        etc, out = _compute_verbose(seq, order)
        return {"ETC": etc, "Trajectory": out}

    # If not verbose, run the compact version and return accordingly
    etc = _compute_compact(seq, order)
    return {"ETC": etc}


def compute_save(seq, filename, order=2):
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
    # Run the verbose version
    etc, out = _compute_verbose(seq, order)

    # Save the output to a csv file and return
    save(out, filename)
    return {"ETC": etc}
