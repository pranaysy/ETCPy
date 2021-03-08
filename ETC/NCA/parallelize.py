#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from ETC.CCMC.pairs import CCM_causality as CCM_compute
from ETC.CCC.compute_CCC import compute as CCC_compute
from multiprocessing import Pool
from functools import partial
from itertools import combinations


def _kernel(inputs, CCC_params):
    """
    Wrapper for computing causality estimates on a sequence pair

    Used for causal discovery and estimation from CCM based methods as well as CCC.

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
   CCC_params : dict
        The following 3 parameters for CCC as key-value pairs:
        "LEN_past" : int
            Parameter "L": Window length of immediate past values of seq_x and seq_y.
        "ADD_meas" : int
            Parameter "w": Window length of present values of seq_x. Minimal data length
            over which CC rate can be reliably estimated, application/domain-specific
        "STEP_size" : int
            Parameter "delta": Step-size for sliding chunks across both sequences. An overlap
            of 20-50% between successive chunks or windows suggested.
        The dictionary can be generated interactively using CCC.get_params()

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

    # Execute the CCM_compute on the sequence pair
    out.update(CCM_compute(seq_x, seq_y))

    # Execute CCC_compute on the sequence pair in one direction
    out.update({"CCC_y_to_x": CCC_compute(seq_x, seq_y, **CCC_params)})

    # Execute CCC_compute on the sequence pair in the other direction
    out.update({"CCC_x_to_y": CCC_compute(seq_y, seq_x, **CCC_params)})

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


def parallelized(pairs, CCC_params):
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

    exec_kernel = partial(_kernel, CCC_params=CCC_params)

    # Initialize pool of parallel workers
    pool = Pool()

    # Confirm to stdout
    print("Running CCM & CCC in parallel on input ... ", end="")

    # Map-execute function across sequences
    out = pool.map_async(exec_kernel, enumerate(pairs))

    # Graceful exit
    pool.close()
    pool.join()

    # Confirm completion
    print("Done!")

    # Return collected results
    return out.get()