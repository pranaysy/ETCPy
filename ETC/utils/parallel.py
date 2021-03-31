#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from joblib import Parallel, delayed
from itertools import islice, combinations, product
import pandas as pd


def get_slices(sequence, chunksize, offset):
    """
    This function takes an input sequence and produces chunks of chosen size.
    Offset can be used to control degree of overlap (or distance between chunks
    that don't overlap)

    Parameters
    ----------
    sequence : tuple or list
        Sequence of integers.
    chunksize : int
        Length of each produced chunk.
    offset : int, optional
        Number of elements to shift each chunk by. The default is 1.
        Setting this to any value less than size allows control of overlap.
        Setting this >= size produces non-overlapping chunks.

    Returns
    -------
    enumerate
        enumerate object that produces indexed chunks of specified size, one at a time.

    """
    if not isinstance(chunksize, int):
        raise TypeError("Argument chunksize should be an int")
    if not isinstance(offset, int):
        raise TypeError("Argument offset should be an int")
    if chunksize >= len(sequence):
        raise ValueError("Argument chunksize must be smaller than sequence length")
    if offset >= len(sequence):
        raise ValueError("Argument offset must be smaller than sequence length")

    return enumerate(
        zip(*(islice(sequence, i, None, offset) for i in range(chunksize)))
    )


def get_crosspairs(collection_a, collection_b):
    """
    Takes two sets of sequences and returns a cross-product of all sequences from both
    sets for bivariate computations. Useful for studying interactions between sets of
    variables.

    Parameters
    ----------
    collection_a : list or generator of sequences
        Multiple sequences belonging to one category or group.
    collection_b : list or generator of sequences
        Multiple sequences belonging to another category or group.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return enumerate(product(collection_a, collection_b))


def get_rowpairs(matrix):
    """
    Create a generator for iterating over pairs of rows of an input matrix

    Parameters
    ----------
    matrix : numpy array, uint32, 2D
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
    # Check for matrix type
    for row1, row2 in combinations(range(0, matrix.shape[0]), 2):
        yield (row1, row2, matrix[row1, :], matrix[row2, :])


def compute(collection, func):
    return pd.DataFrame(
        Parallel(n_jobs=-1, verbose=9)(delayed(func)(record) for record in collection)
    )
