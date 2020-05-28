#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter
from math import log2
from random import choices
from ETC.seq import estimates, recode
from array import array
import re

def sanitize(text, whitespace=False, lowercase=False):

    if whitespace:
        joiner = " "

    else:
        joiner = ""

    if lowercase:
        text = text.lower()

    text = joiner.join(re.findall("[a-zA-Z]+", text))

    return text


# Function definitions
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
    # Get smallest value
    a = min(seq)

    # Compute reciprocal of peak-to-peak per bin
    delta_inv = n_bins / (max(seq) - a + 1e-6)

    # Transform each element and return
    return [1 + int((elem - a) * delta_inv) for elem in seq]


def generate(size=10, partitions=2):
    """
    This function generates discrete random data of desired size and bins.

    Parameters
    ----------
    size : int, optional
        Length of sequence to generate. The default is 10.
    partitions : int, optional
        Number of bins/paritions to create.

    Returns
    -------
    list
        Collection of integers sampled from discrete uniform.

    """
    if not (isinstance(partitions, int) and isinstance(size, int) and partitions >= 2):
        print(">> Number of bins is invalid ...")
        return None

    return recode.cast(choices(range(1, partitions + 1), k=size))

def frequencies(seq):

    return Counter(seq).most_common()

def entropy(seq, legacy=False):
    """
    This function computes Shannon Entropy of a given sequence.

    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.

    Returns
    -------
    float
        Shannon entropy of sequence.

    """

    if isinstance(seq, array) and seq.typecode == 'I' and not legacy:
        return estimates.entropy(seq)

    # Get counts from Counter, normalize by total, transform each and sum all
    return sum(
        -seq * log2(seq) for seq in (elem / len(seq) for elem in Counter(seq).values())
    )

