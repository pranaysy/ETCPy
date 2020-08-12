#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""
from collections import Counter
from math import log2
from random import choices
from random import seed as seedvalue
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


def generate(size=10, partitions=2, seed=None):
    """
    This function generates discrete random data of desired size and bins.

    Parameters
    ----------
    size : int, optional
        Length of sequence to generate. The default is 10.
    partitions : int, optional
        Number of bins/paritions to create.
    seed : int, optional
        Seed value for initializing the random number generator. The default is None

    Returns
    -------
    list
        Collection of integers sampled from discrete uniform.

    """
    if not (isinstance(partitions, int) and isinstance(size, int) and partitions >= 2):
        print(partitions, size)
        print(">> Number of bins is invalid ...")
        return None

    if seed:
        seedvalue(seed)

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

    if isinstance(seq, array) and seq.typecode == "I" and not legacy:
        return estimates.entropy(seq)

    # Get counts from Counter, normalize by total, transform each and sum all
    return sum(
        -seq * log2(seq) for seq in (elem / len(seq) for elem in Counter(seq).values())
    )
