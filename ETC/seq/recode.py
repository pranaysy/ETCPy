#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from string import ascii_lowercase
from random import shuffle, choices
from array import array
import numpy as np
from ETC.seq.check import zeroes


def cast(seq):

    if seq is not None:
        if isinstance(seq, np.ndarray):
            try:
                out = seq.astype("uint32")
                if zeroes(out):
                    print("> Input contains 0!")
                    print("> Symbols shifted up by 1 ")
                    return out + 1
                return out

            except TypeError as error:
                print("ERROR:", error)
                print("> Input must be a list/tuple/array of positive integers!")
                print('> Recode or partition using "ETC.seq.recode"')
                return None

            except OverflowError as error:
                print("ERROR:", error)
                print("> Input must be a list/tuple/array of positive integers!")
                print('> Recode or partition using "ETC.seq.recode"')
                return None

        else:
            try:
                out = array("I", seq)
                if zeroes(out):
                    print("> Input contains 0!")
                    print('> Recode or partition using "ETC.seq.recode" ')
                    return None
                return out

            except TypeError as error:
                print("ERROR:", error)
                print("> Input must be a list/tuple/array of positive integers!")
                print('> Recode or partition using "ETC.seq.recode"')
                return None

            except OverflowError as error:
                print("ERROR:", error)
                print("> Input must be a list/tuple/array of positive integers!")
                print('> Recode or partition using "ETC.seq.recode"')
                return None

    print("No input sequence provided.")
    return None


def recode_lexical(text, case_sensitive=True):

    if not isinstance(text, str):
        print("ERROR: Input is not a string.")
        return None
    if not case_sensitive:
        text = text.lower()
    alphabets = sorted(set(text))
    replacer = dict((y, x + 1) for x, y in enumerate(alphabets))
    text = cast([replacer[x] for x in text])
    return text


def recode_alphabetical(text):

    text = text.lower()
    if not set(text).issubset(ascii_lowercase):
        print("> Input contains non alphabetical characters!")
        return None
    replacer = dict((y, x + 1) for x, y in enumerate(ascii_lowercase))
    text = cast([replacer[x] for x in text])
    return text


def recode_dna(text):

    replacer = {"A": 1, "G": 1, "C": 2, "T": 2}
    text = cast([replacer[x] for x in text.upper()])
    return text


def recode_random(text):

    alphabets = list(set(text))
    shuffle(alphabets)
    replacer = dict((y, x + 1) for x, y in enumerate(alphabets))
    text = cast([replacer[x] for x in text])
    return text


def recode_randint(text):

    alphabets = list(set(text))
    numbers = choices(range(1, 2 ** 20), k=len(alphabets))
    replacer = dict(zip(alphabets, numbers))
    text = cast([replacer[x] for x in text])
    return text


def partition(seq, n_bins):
    """
    This function takes an input sequence and bins it into discrete points.

    Parameters
    ----------
    seq : list/tuple of float
        Collection of floats.
    n_bins : int
        Number of bins/paritions to create.

    Returns
    -------
    list
        Collection of integers. Contains unique integers from 1 to n_bins.

    """
    assert (
        isinstance(n_bins, int) and n_bins > 1
    ), "ERROR: Number of bins should be a positive integer"

    # Get smallest value
    a = min(seq)

    # Compute reciprocal of peak-to-peak per bin
    delta_inv = n_bins / (max(seq) - a + 1e-6)

    # Transform each element and return
    return [1 + int((elem - a) * delta_inv) for elem in seq]


def partition_numpy(nparr, n_bins):
    """
    This function takes an input sequence & partitions it into equiwidth discrete bins.

    Min-max scaling, followed by equiwidth binning for each row

    Parameters
    ----------
    nparr : numpy array, int or float, 2D
        Each row representing a different sequence. (Columns as time)
    n_bins : int
        Number of bins/paritions to create.

    Returns
    -------
    list
        Collection of integers. Contains unique integers from 1 to n_bins.

    """
    assert (
        isinstance(n_bins, int) and n_bins > 1
    ), "ERROR: Number of bins should be a positive integer"

    assert (
        isinstance(nparr, np.ndarray) and nparr.ndim == 2
    ), ">ERROR: Input must be 2D NumPy array of numbers"

    # Get smallest value
    a = nparr.min(axis=1)[:, np.newaxis]

    # Compute reciprocal of peak-to-peak per bin
    delta_inv = n_bins / (nparr.max(axis=1)[:, np.newaxis] - a + 1e-6)

    # Transform each element and return
    return 1 + ((nparr - a) * delta_inv).astype("uint32")
