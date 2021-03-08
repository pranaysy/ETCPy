#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the Compression-Complexity based Causality for two sequences.

@author: Pranay S. Yadav
"""

# Import functions
from functools import partial
from ETC import compute_1D, compute_2D
from ETC.seq.recode import partition, cast
from ETC.seq.check import arraytype

import numpy as np
import array

# Curry the functions for computing 1D and 2D ETC estimates
get1D = partial(compute_1D, order=2, verbose=False, truncate=True)
get2D = partial(compute_2D, order=2, verbose=False, truncate=True)


def get_params():
    """
    Helper function for creating a dictionary of CCC params interactively

    Returns
    -------
    dict
        The main 3 parameters for passing directly into CCC.compute

    """
    print("#" * 80)
    print("-" * 80)
    print("Initialize CCC Parameters: The thorns")
    print("-" * 80, "\n")
    print("All the following 3 should be integers\n")
    LEN_past = int(input("1. Window length of immediate past values: "))
    ADD_meas = int(input("2. Window length of present values: "))
    STEP_size = int(input("3. Step-size for sliding window ahead: "))
    print("\n", "-" * 80)
    print("#" * 80)
    return {"LEN_past": LEN_past, "ADD_meas": ADD_meas, "STEP_size": STEP_size}


def compute(seq_x, seq_y, LEN_past, ADD_meas, STEP_size, n_partitions=False):
    """
    Estimate the Compression-Complexity based Causality for two sequences.

    The direction of causality being assessed is from seq_y -> seq_x. Various other
    parameters need to be specified, a brief description is offered below.

    For detailed explanations regarding the parameters, interpretations and of the inner
    workings, please refer to the research article along with the supplementary:
        Kathpalia, Aditi, and Nithin Nagaraj. “Data-Based Intervention Approach for
        Complexity-Causality Measure.” PeerJ Computer Science 5 (May 2019): e196.
        https://doi.org/10.7717/peerj-cs.196.

    Parameters
    ----------
    seq_x : list or tuple
        Sequence of numbers, if not integers specify n_partitions for binning.
    seq_y : list or tuple
        Sequence of numbers, if not integers specify n_partitions for binnings.
    LEN_past : int
        Parameter "L": Window length of immediate past values of seq_x and seq_y.
    ADD_meas : int
        Parameter "w": Window length of present values of seq_x. Minimal data length
        over which CC rate can be reliably estimated, application/domain-specific
    STEP_size : int
        Parameter "delta": Step-size for sliding chunks across both sequences. An overlap
        of 20-50% between successive chunks or windows suggested.
    n_partitions : int or bool, optional
        Parameter "B": Number of bins. Smalles number of symbols that capture the time
        series dynamics. The default is False indicating that the data is already in the
        form of discrete symbolic sequences.

    Returns
    -------
    CCC : float
        Estimated Compression-Complexity based Causality for direction seq_y -> seq_x.

    """
    # Sanity checks
    assert len(seq_x) == len(seq_y), "ERROR: Sequences must have the same length!"
    assert (
        isinstance(LEN_past, int) and LEN_past > 1
    ), "ERROR: LEN_past must be a positive integer!"
    assert (
        isinstance(ADD_meas, int) and ADD_meas > 1
    ), "ERROR: ADD_meas must be a positive integer!"
    assert (
        isinstance(STEP_size, int) and STEP_size > 1
    ), "ERROR: STEP_size must be a positive integer!"

    # Partition data if requested with the specificed number of bins
    if n_partitions:
        seq_x = partition(seq_x, n_partitions)
        seq_y = partition(seq_y, n_partitions)

    # Check whether input is a discrete symbolic sequence
    if not arraytype(seq_x):
        seq_x = cast(seq_x)
    if not arraytype(seq_y):
        seq_y = cast(seq_y)

    # Set switch for operating differently on native vs numpy arrays
    if type(seq_x) == np.ndarray or type(seq_y) == np.ndarray:
        combine = lambda x, y: np.hstack([x, y])
    if type(seq_x) == array.array or type(seq_y) == array.array:
        combine = lambda x, y: x + y

    # Setup variables
    LEN = len(seq_x)
    LEN_to_check = LEN_past + ADD_meas

    # Initialize aggregators
    l_1D = []
    l_2D = []

    # Iterate over chunks of both sequences
    for k in range(0, LEN - LEN_to_check, STEP_size):

        ## Compression-Complexity of past values of seq_x
        # 1D ETC of a chunk of seq_x of length LEN_past
        ETC1D_ini = get1D(seq_x[k : k + LEN_past])["NETC1D"]

        ## Compression-Complexity of past values of seq_x and seq_y
        # 2D ETC of chunks of both seq_x,seq_y of length LEN_past at the same locus
        ETC2D_ini = get2D(seq_x[k : k + LEN_past], seq_y[k : k + LEN_past],)["NETC2D"]

        ## Compression-Complexity of present values of seq_x
        # 1D ETC of a chunk of seq_x of length LEN_to_check
        ETC1D_fin = get1D(seq_x[k : k + LEN_to_check])["NETC1D"]

        ## Compression-Complexity of values of seq_x & past of seq_y + present of seq_x
        # 2D ETC of chunks of both seq_x, seq_y of length LEN_to_check at the same locus
        ETC2D_fin = get2D(
            seq_x[k : k + LEN_to_check],
            combine(seq_y[k : k + LEN_past], seq_x[k + LEN_past : k + LEN_to_check]),
        )["NETC2D"]

        # Dynamic Compression-Complexity of seq_x
        ETC1D_delta = ETC1D_fin - ETC1D_ini

        # Dynamic Compression Complexity of seq_x conditional on seq_y
        ETC2D_delta = ETC2D_fin - ETC2D_ini

        # Aggregate Dynamic CCs
        l_1D.append(ETC1D_delta)
        l_2D.append(ETC2D_delta)

    ## Compute Compession-Complexity Causality
    # Average of the difference: CC(X | X_past) - CC(X | Y_past + X_present)
    CCC = (sum(l_1D) - sum(l_2D)) / len(l_1D)
    # print(f"CCC for seq_y -> seq_x = {CCC}")
    return CCC
