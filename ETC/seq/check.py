#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: Pranay S. Yadav
"""

from ETC.seq import estimates
from array import array
import numpy as np


def zeroes(seq):
    if 0 in seq:
        return True
    return False


def equality(seq, legacy=False):
    """
    This function checks if all elements of a collection are equal.
    Parameters
    ----------
    seq : list or tuple
        Sequence of integers.
    Returns
    -------
    bool
        True if all elements equal.
    """

    if arraytype(seq) and not legacy:
        return estimates.equality(seq)

    # Iterate over all elements in sequence
    for element in seq:

        # Break at first inequality
        if seq[0] != element:
            return False

    # Else all equal
    return True


def arraytype(seq):
    if isinstance(seq, array) and seq.typecode == "I":
        return True
    if isinstance(seq, np.ndarray) and seq.dtype == "uint32":
        return True
    return False
