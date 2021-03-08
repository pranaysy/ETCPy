#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from ETC.CCMC.pairs import ETC_causality, LZ_causality, CCM_causality
from multiprocessing import Pool
from functools import partial
from itertools import combinations


def _kernel_seq(inputs, estimator):
    """
    Wrapper around a function that computes anything on two sequences and returns a dict

    While it is written as a general purpose kernel for anything, here it is used for
    causal discovery and estimation from CCM based methods.

    The function unpacks inputs into an index element and a sequence pair and runs the
    estimator function on the sequence pair, returning various estimates in a dict

    Parameters
    ----------
    inputs : tuple
        Tuple of two elements - (a, b) where a is an index, b is a tuple of two. a can
        be produced manually or more typically using enumerate; b holds the two sequences
        usually passed in by zip-ping larger iterables or itertools' product/combinations.
        a, the index, is passed to keep track of order in case of asynchronous execution
        Should look like this: (index, (sequence_x, sequence_y)
    estimator : function
        A function that can compute something on two arrays and return a dict. Preferably
        one that can compute something meaningful, like causal discovery

    Returns
    -------
    out : dict
        Estimates obtained by running estimator on inputs.

    """
    # Unpack inputs
    idx, seqs = inputs

    # Unpack sequences
    idx_x, idx_y, seq_x, seq_y = seqs

    # Initialize dictionary of output estimates with index
    out = {"index_pair": idx, "index_x": idx_x, "index_y": idx_y}

    # Execute the estimator function on the sequence pair
    out.update(estimator(seq_x, seq_y))

    # Some feedback to console
    # print(".", end="")

    return out


def get_rowpairs(matrix):
    """
    Create a generator for iterating over pairs of rows of an input matrix

    Parameters
    ----------
    matrix : numpy array, int or float, 2D
        Each row representing a different sequence. (Columns as time)

    Yields
    ------
    row1 : int
        Index of first row in the pair.
    row2 : int
        Index of second row in the pair.
    np.array, 1D, int
        Data of first row in the pair.
    np.array, 1D, int
        Data of first row in the pair.

    """
    for row1, row2 in combinations(range(0, matrix.shape[0]), 2):
        yield (row1, row2, matrix[row1, :], matrix[row2, :])


def parallelized(pairs, kernel="CCM"):
    """
    This function operates concurrently on a collection of sequence pairs and computes
    estimates using the chosen kernel function.

    Here used for computing causal estimates from sequences pairs in batch, each pair
    runs on a separate CPU core as a process.

    CAUTION: main module is unguarded, do not run these functions as is,
        particularly on Windows!

    Parameters
    ----------
    pairs : list/tuple/generator
        Collection of pairs of integer sequences.
    kernel : str, optional
        Name of an estimator function. Currently available: "CCM", "ETC" and "LZ". The
        default is "CCM".

    Returns
    -------
    list of dict elements
        Each dictionary element contains index, length of sequence & ETC.

    """

    if kernel == "CCM":
        exec_kernel = partial(_kernel_seq, estimator=CCM_causality)
    elif kernel == "ETC":
        exec_kernel = partial(_kernel_seq, estimator=ETC_causality)
    elif kernel == "LZ":
        exec_kernel = partial(_kernel_seq, estimator=LZ_causality)
    else:
        print("> ERROR: Invalid kernel specified")
        return None

    # Initialize pool of parallel workers
    pool = Pool()

    # Confirm to stdout
    print(f"Running kernel={kernel} in parallel on input ... ", end="")

    # Map-execute function across sequences
    out = pool.map_async(exec_kernel, enumerate(pairs))

    # Graceful exit
    pool.close()
    pool.join()

    # Confirm completion
    print("Done!")

    # Return collected results
    return out.get()
