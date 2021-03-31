#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for computing markov transition probability matrices from
nucleotide sequence stored as text files.

compute() is the main function that wraps around smaller modular functions.

@author: pranay
"""

# Import calls
from pathlib import Path
from random import choices, seed
import numpy as np
import pandas as pd

# Function Definitions
def _read_sequence(filepath):
    """
    This function reads a file & returns it as a string.
    Uses pathlib's functionality & is called by the wrapper compute()

    Parameters
    ----------
    filepath : Path object
        Valid path to file containing nucleotide sequence.

    Returns
    -------
    string
        String containing nucleotide sequence.

    """
    # Use Path object's hook to read file as text and return
    return filepath.read_text()


def _generate_overlaps(sequence, order):
    """
    This function takes an input sequence & generates overlapping subsequences
    of length order + 1, returned as a tuple. Has no dependencies & is called
    by the wrapper compute()

    Parameters
    ----------
    sequence : string
        String containing nucleotide sequence.
    order : int
        Order of Markov Transition Probability Matrix for computing overlaps

    Returns
    -------
    tuple
        Contains all overlapping sub-sequences (of length order + 1) for given
        order. Total number of sub-sequences is len(sequence) - (order + 1)

    """
    # Initialize aggregator
    aggregator = []

    # Increment order by 1 as the Markov model includes the current state, such
    # that the length of sequence corresponding to a state is order + 1
    order += 1

    # Iteratively store sequences shifted to the left by 1 step
    for idx in range(order):
        aggregator.append(sequence[idx : idx - order])

    # Join the shifted sequences through element-wise concatenation & return
    return tuple(map("".join, zip(*aggregator)))


def _compute_transition_probs(sequence, compact=True, flatten=False):
    """
    This function takes a tuple of strings (overlapping subsequences) as input
    and computes transition probabilities, returned as a dataframe. Switches
    control the form (compactness & wide-form vs long-form) of output dataframe.

    Parameters
    ----------
    sequence : string
        String containing nucleotide sequence.
    compact : bool, optional
        Whether to return the full sparse matrix or to return a more compact
        representation of it. If False, returns the a square matrix with most
        elements zero. If True, returns non-zero columns only for all rows.
            The default is True.
    flatten : bool, optional
        Whether to flatten (or tidy / long-form) the matrix or not (wide-form).
        If True, returns only one column containing probabilities through a
        row-wise representation. If False, returns multiple columns with
        probabilities.
            The default is False.

    Returns
    -------
    pandas DataFrame
        Tabulated transition probabilities with row & column labels describing
        the (n-1)th and nth state respectively. If flatten is True, column
        labels for nth state are transposed into a column such that there are
        2 columns, 1 each for the (n-1)th and nth state.

    """
    # If compact requested, use only the last alphabet of the next subsequence.
    # The Nth element will only differ from the (N-1)th in shift by 1
    if compact:

        # Convert tuple to numpy array, excluding last element
        temp = np.array(sequence[:-1])

        # Extract the last alphabet from each word except the first
        next_alphabet = np.array([x[-1] for x in sequence[1:]])

        # Compute normalized frequencies via cross-tabulation
        df = pd.crosstab(temp, next_alphabet, normalize="index")

    # If full matrix to be returned, cross-tabulate shifted sequences
    else:

        # Convert tuple containing overlapping subsequences to numpy array
        temp = np.array(sequence)

        # Compute normalized frequencies via cross-tabulation
        df = pd.crosstab(temp[:-1], temp[1:], normalize="index")

    # Set proper identifier labels
    df.index.name = "previous"
    df.columns.name = "next"

    # If flatten requested, pivot all columns into a single column
    if flatten:

        # Stack all columns
        df = df.stack()

        # Set name for the Series object
        df.name = "probability"

        # Return a DataFrame by resetting the Series index
        return df.reset_index()

    # If flatten not requested, return DataFrame
    return df


def _check_inputs(filepath, order, compact, flatten):
    """
    This function checks the input arguments to compute() for validity based
    on descriptions below.

    Parameters
    ----------
    filepath : Path object
        Valid path to file containing nucleotide sequence.
    order : int
        Order of Markov Transition Probability Matrix for computing overlaps
    compact : bool, optional
        Whether to return the full sparse matrix or to return a more compact
        representation of it
    flatten : bool, optional
        Whether to flatten (or tidy / long-form) the matrix or not.

    Returns
    -------
    bool
        True if all inputs are valid.

    """
    # Check type of input path
    if not isinstance(filepath, Path):
        print("> ERROR: Input should be a Path object ...")
        return False

    # Check if path exists and points to a file
    if not (filepath.exists() and filepath.is_file()):
        print("> ERROR: Path does not exist ...")
        return False

    # Check if order is a non-negative integer
    if not (isinstance(order, int) and order >= 0):
        print("> ERROR: order should be a non-negative integer ...")
        return False

    # Check if other args are boolean types
    if not (isinstance(compact, bool) and isinstance(flatten, bool)):
        print("> ERROR: compact and flatten args should be a boolean ...")
        return False

    # If all inputs are valid, yay
    return True


def compute(filepath, order, compact=True, flatten=False):
    """
    This function takes a text file containing a nucleotide sequence, computes
    computes the transition probability matrix of given order and returns it as
    a labelled dataframe. Output can be tuned through optional switches. This
    function is modular and wraps around the following 4 functions:
        _check_inputs - for validating input arguments
        _read_sequence - for reading text file containing nucleotide sequence
        _generate_overlaps - for creating overlapped subsequences
        _compute_transition_probs - for creating transition probability matrix

    Parameters
    ----------
    filepath : Path object
        Valid path to file containing nucleotide sequence.
    order : int
        Order of Markov Transition Probability Matrix for computing overlaps
    compact : bool, optional
        Whether to return the full sparse matrix or to return a more compact
        representation of it. If False, returns the a square matrix with most
        elements zero. If True, returns non-zero columns only for all rows.
            The default is True.
    flatten : bool, optional
        Whether to flatten (or tidy / long-form) the matrix or not. If True,
        returns only one column containing probabilities through a row-wise
        representation. If False, returns multiple columns with probabilities.
            The default is False.

    Returns
    -------
    pandas DataFrame
        Tabulated transition probabilities with row & column labels describing
        the (n-1)th and nth state respectively. If flatten is True, column
        labels for nth state are transposed into a column such that there are
        2 columns, 1 each for the (n-1)th and nth state.

    """
    # If any input is not valid, break
    if not _check_inputs(filepath, order, compact, flatten):
        return None

    # Read sequence file and get tuple of overlapping subsequences
    sequence = _generate_overlaps(_read_sequence(filepath), order)

    # Compute and return transition probability matrix
    return _compute_transition_probs(sequence, compact=compact, flatten=flatten)


def sample_sequence(sequence, order, size, sampler_seed=0):

    # Read sequence file and get tuple of overlapping subsequences
    overlapped = _generate_overlaps(sequence, order)

    # Compute and return transition probability matrix
    transition_probs = _compute_transition_probs(
        overlapped, compact=True, flatten=False
    )

    order += 1
    chain = sequence[-order:]

    seed(sampler_seed)

    for n in range(size):
        last = chain[-order:]
        probs = transition_probs.loc[last, :]
        new = "".join(choices(population=probs.index, weights=probs.values, k=1))
        chain += new

    return chain[order:]
